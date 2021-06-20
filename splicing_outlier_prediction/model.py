from tqdm import tqdm
import pandas as pd
try:
    from mmsplice import MMSplice
    from mmsplice.utils import df_batch_writer, delta_logit_PSI_to_delta_PSI
except ImportError:
    pass
from splicing_outlier_prediction.result import SplicingOutlierResult


class SpliceOutlier:

    def __init__(self, clip_threshold=None):
        import mmsplice
        self.mmsplice = MMSplice()
        self.clip_threshold = clip_threshold

    def _add_delta_psi_single_ref(self, df, splicemap):
        df_splicemap = splicemap.df

        core_cols = ['junctions', 'Chromosome', 'Start', 'End', 'Strand']
        cols_splicemap = df_splicemap.columns.difference(core_cols, False)
        df_splicemap = df_splicemap.rename(columns={'junctions': 'junction'}) \
            .set_index('junction')[cols_splicemap]

        df = df[df.columns.difference(cols_splicemap, False)] \
            .set_index('junction')

        df_joined = df.join(df_splicemap, how='inner').reset_index()

        delta_psi = delta_logit_PSI_to_delta_PSI(
            df_joined['delta_logit_psi'],
            df_joined['ref_psi'],
            clip_threshold=self.clip_threshold or 0.01
        )
        df_joined.insert(8, 'delta_psi', delta_psi)
        df_joined.insert(3, 'tissue', splicemap.name)
        return df_joined

    def _add_delta_event(self, df, dl, event_type):
        df = df[df['event_type'] == event_type]
        return pd.concat(
            self._add_delta_psi_single_ref(df, splicemap)
            for splicemap in dl.splicemaps5
        )

    def _add_delta_psi(self, df, dl):
        return pd.concat([
            self._add_delta_event(df, dl, 'psi5'),
            self._add_delta_event(df, dl, 'psi3')
        ])

    def predict_on_batch(self, batch, dataloader):
        columns = batch['metadata']['junction'].keys()
        df = self.mmsplice._predict_batch(batch, columns)
        del df['exons']
        df = df.rename(columns={'ID': 'variant'})
        df_with_delta_psi = self._add_delta_psi(df, dataloader)
        return df_with_delta_psi

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

    def predict_save(self, dataloader, output_csv,
                     batch_size=512, progress=True):
        df_batch_writer(self._predict_on_dataloader(dataloader), output_csv)
