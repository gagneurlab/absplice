import itertools
from tqdm import tqdm
import pandas as pd
from kipoi.data import SampleIterator
try:
    from mmsplice import MMSplice
    from mmsplice.junction_dataloader import JunctionPSI5VCFDataloader, \
        JunctionPSI3VCFDataloader
    from mmsplice.utils import encodeDNA, df_batch_writer, \
        delta_logit_PSI_to_delta_PSI
except ImportError:
    pass
from splicing_outlier_prediction.dataloader_updated import RefTableMixin
from splicing_outlier_prediction.result import SplicingOutlierResult


class SpliceOutlierDataloader(RefTableMixin, SampleIterator):

    def __init__(self, fasta_file, vcf_file, ref_tables5=list(), ref_tables3=list(), **kwargs):
        super().__init__(ref_tables5, ref_tables3, **kwargs)
        import mmsplice
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self._generator = iter([])

        if self.intron_annotation5:
            self.dl5 = JunctionPSI5VCFDataloader(
                self.intron_annotation5, fasta_file, vcf_file, encode=False, **kwargs)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl5, self.intron_annotation5, event_type='psi5'))

        if self.ref_tables3:
            self.dl3 = JunctionPSI3VCFDataloader(
                self.intron_annotation3, fasta_file, vcf_file, encode=False, **kwargs)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl3, self.intron_annotation3, event_type='psi3'))

    def _iter_dl(self, dl, intron_annotations, event_type):
        for row in dl:
            junction_id = row['metadata']['exon']['junction']
            ref_row = intron_annotations.df.loc[junction_id]
            row['metadata']['junction'] = dict()
            row['metadata']['junction']['junction'] = ref_row.name
            row['metadata']['junction']['event_type'] = event_type
            row['metadata']['junction'].update(ref_row.to_dict())
            yield row

    def __next__(self):
        return next(self._generator)

    def __iter__(self):
        return self

    def batch_iter(self, batch_size=32, **kwargs):
        for batch in super().batch_iter(batch_size, **kwargs):
            batch['inputs']['seq'] = self._encode_batch_seq(
                batch['inputs']['seq'])
            batch['inputs']['mut_seq'] = self._encode_batch_seq(
                batch['inputs']['mut_seq'])
            yield batch

    def _encode_batch_seq(self, batch):
        return {k: encodeDNA(v.tolist()) for k, v in batch.items()}


class SpliceOutlier:

    def __init__(self, clip_threshold=None):
        import mmsplice
        self.mmsplice = MMSplice()
        self.clip_threshold = clip_threshold

    def _add_delta_psi(self, df):
        df_list = list()
        columns = ['junctions', 'Chromosome', 'Start', 'End', 'Strand']
        # Go through all ref_tables and add delta_psi
        for ref_table in self.ref_tables:
            # TODO: get tissue information and store into column
            col = [c for c in set(ref_table.columns).difference(set(columns))]
            df_joined = ref_table.set_index('junctions')[col].join(df.set_index('junctions'))
            delta_psi = delta_logit_PSI_to_delta_PSI(
                df_joined['delta_logit_psi'],
                df_joined['ref_psi'],
                clip_threshold=self.clip_threshold or 0.01
            )
            df_joined.insert(18, 'delta_psi', delta_psi)
            df_list.append(df_joined)
        df_with_delta_psi = pd.concat(df_list)
        return df_with_delta_psi

    def predict_on_batch(self, batch):
        columns = [
            'samples', 'maf', 'genotype', 'GQ', 'DP_ALT',
            *batch['metadata']['junction'].keys()
        ]
        df = self.mmsplice._predict_batch(batch, columns)
        del df['exons']
        df = df.rename(columns={'ID': 'variant'})
        df_with_delta_psi = self._add_delta_psi(df)
        return df_with_delta_psi

    def _predict_on_dataloader(self, dataloader,
                               batch_size=512, progress=True):
        # TODO: QUESTION: how does it know that this dataloader comes from class above?
        dt_iter = dataloader.batch_iter(batch_size=batch_size)
        if progress:
            dt_iter = tqdm(dt_iter)

        for batch in dt_iter:
            yield self.predict_on_batch(batch)

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
