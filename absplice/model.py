from tqdm import tqdm
import pandas as pd
try:
    from mmsplice import MMSplice
    from mmsplice.utils import df_batch_writer, df_batch_writer_parquet, delta_logit_PSI_to_delta_PSI
except ImportError:
    pass
from absplice.result import SplicingOutlierResult
from pathlib import Path
import pathlib


class SpliceOutlier:

    def __init__(self, clip_threshold=None):
        import mmsplice
        self.mmsplice = MMSplice()
        self.clip_threshold = clip_threshold

    def _add_metadata_event(self, df, metadata, event_type):
        df = df[df['event_type'] == event_type].set_index('junction')
        
        return df.join(pd.DataFrame(
            [
                row
                for junc in df.index
                for row in metadata[junc]
            ], 
            columns=['junction', 'gene_id', 'tissue', 'ref_psi', 'median_n', 'gene_name', 'splice_site']
        ).set_index('junction')).reset_index().drop_duplicates()

    def _add_metadata(self, df, dl):
        return pd.concat([
            self._add_metadata_event(df, dl.metadata_splicemap5, 'psi5'),
            self._add_metadata_event(df, dl.metadata_splicemap3, 'psi3')
        ])
        
    def _add_delta_psi(self, df):
        delta_psi = delta_logit_PSI_to_delta_PSI(
            df['delta_logit_psi'],
            df['ref_psi'],
            clip_threshold=self.clip_threshold or 0.01
        )
        df.insert(8, 'delta_psi', delta_psi)
        return df

    def predict_on_batch(self, batch, dataloader):
        columns = batch['metadata']['junction'].keys()
        df = self.mmsplice._predict_batch(batch, columns)
        del df['exons']
        df = df.rename(columns={'ID': 'variant'})
        df = self._add_metadata(df, dataloader)
        df = self._add_delta_psi(df)
        cols = [
            'variant', 'tissue', 'junction', 'event_type',
            'splice_site', 'ref_psi', 'median_n', 
            'gene_id', 'gene_name',
            'delta_logit_psi', 'delta_psi',
        ]
        df = df[cols]
        return df

    def _predict_on_dataloader(self, dataloader,
                               batch_size=512, progress=True):
        dt_iter = dataloader.batch_iter(batch_size=batch_size)
        if progress:
            dt_iter = tqdm(dt_iter)

        for batch in dt_iter:
            yield self.predict_on_batch(batch, dataloader)

    def predict_on_dataloader(self, dataloader, batch_size=512, progress=True):
        return SplicingOutlierResult(pd.concat(
            self._predict_on_dataloader(
                dataloader,
                batch_size=batch_size,
                progress=progress)
        ))

    def predict_save(self, dataloader, output_path,
                     batch_size=512, progress=True):
        if not isinstance(output_path, pathlib.PosixPath):
            output_path = Path(output_path)
        if output_path.suffix.lower() == '.csv':
            df_batch_writer(self._predict_on_dataloader(dataloader), output_path)
        elif output_path.suffix.lower() == '.parquet':
            df_batch_writer_parquet(self._predict_on_dataloader(dataloader), output_path)
