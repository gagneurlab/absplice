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

    def _add_delta_psi_single_ref(self, df, ref_table):
        cols_ref_table = list(set(ref_table.columns).difference(
            set(['junctions', 'Chromosome', 'Start', 'End', 'Strand'])))
        ref_table = ref_table.rename(columns={'junctions': 'junction'})
        cols_df = [c for c in df.columns if (c not in cols_ref_table)]
        df_joined = ref_table.set_index('junction')[cols_ref_table].join(
            df[cols_df].set_index('junction'), how='inner').reset_index()
        delta_psi = delta_logit_PSI_to_delta_PSI(
            df_joined['delta_logit_psi'],
            df_joined['ref_psi'],
            clip_threshold=self.clip_threshold or 0.01
        )
        df_joined.insert(18, 'delta_psi', delta_psi)
        return df_joined

    def _add_delta_psi5(self, df, dataloader):
        df5_list = list()
        df5 = df[df['event_type'] == 'psi5']
        for ref_table5 in dataloader.combined_ref_tables5.ref_tables:
            df5_joined = self._add_delta_psi_single_ref(df5, ref_table5)
            df5_list.append(df5_joined)
        df5_with_delta_psi = pd.concat(df5_list)
        return df5_with_delta_psi

    def _add_delta_psi3(self, df, dataloader):
        df3_list = list()
        df3 = df[df['event_type'] == 'psi3']
        for ref_table3 in dataloader.combined_ref_tables3.ref_tables:
            df3_joined = self._add_delta_psi_single_ref(df3, ref_table3)
            df3_list.append(df3_joined)
        df3_with_delta_psi = pd.concat(df3_list)
        return df3_with_delta_psi

    def _add_delta_psi(self, df, dataloader):
        df5_with_delta_psi = self._add_delta_psi5(df, dataloader)
        df3_with_delta_psi = self._add_delta_psi3(df, dataloader)
        df_with_delta_psi = pd.concat([df5_with_delta_psi, df3_with_delta_psi])
        return df_with_delta_psi

    def predict_on_batch(self, batch, dataloader):
        columns = [
            'samples', 'maf', 'genotype', 'GQ', 'DP_ALT',
            *batch['metadata']['junction'].keys()
        ]
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
