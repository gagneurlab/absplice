from collections import defaultdict
import itertools
from typing import List
import pandas as pd
from tqdm import tqdm
from kipoi.data import SampleIterator
from splicemap.splice_map import SpliceMap
from absplice.utils import junction_str_to_tuple

try:
    from mmsplice.junction_dataloader import JunctionPSI5VCFDataloader, \
        JunctionPSI3VCFDataloader
    from mmsplice.utils import encodeDNA
except ImportError:
    pass


class SpliceMapMixin:

    def __init__(self, splicemap5=None, splicemap3=None, progress=True):
        self.progress = progress

        if splicemap5 is None and splicemap3 is None:
            raise ValueError(
                '`ref_tables5` and `ref_tables3` cannot be both empty')

        if splicemap5 is not None:
            self.splicemaps5 = self._read_splicemap(splicemap5)
            self.metadata_splicemap5 = self._splicemap_metadata(self.splicemaps5)
            self.combined_splicemap5 = self.metadata_splicemap5.index.unique()
        else:
            self.combined_splicemap5 = None

        if splicemap3 is not None:
            self.splicemaps3 = self._read_splicemap(splicemap3)
            self.metadata_splicemap3 = self._splicemap_metadata(self.splicemaps3)
            self.combined_splicemap3 = self.metadata_splicemap3.index.unique()
        else:
            self.combined_splicemap3 = None

    def _splicemap_metadata(self, splicemaps):
        dfs = [
            (
                sm.df
                .astype({
                    "Chromosome": pd.StringDtype(),
                    "Start": pd.Int32Dtype(),
                    "End": pd.Int32Dtype(),
                    "Strand": pd.StringDtype()}
                )
                .assign(tissue=sm.name)
                .set_index(["Chromosome", "Start", "End", "Strand"])
            ) for sm in splicemaps
        ]
        if len(splicemaps) > 1:
            df = pd.concat(dfs)
        else:
            df = dfs[0]

        return df

    @staticmethod
    def _read_splicemap(path):
        if type(path) is str:
            return [SpliceMap.read_csv(path)]
        elif type(path) is SpliceMap:
            return [path]
        elif type(path) is list:
            return [SpliceMapMixin._read_splicemap(i)[0] for i in path]
        else:
            print(type(path))
            raise ValueError(
                '`splicemap5` or `splicemap3` arguments should'
                ' be list of path to splicemap files'
                ' or `SpliceMap` object')


class SpliceOutlierDataloader(SpliceMapMixin, SampleIterator):

    def __init__(self, fasta_file, vcf_file, splicemap5=None, splicemap3=None):
        SpliceMapMixin.__init__(self, splicemap5, splicemap3)

        import mmsplice
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self._generator = iter([])

        if self.combined_splicemap5 is not None:
            self.dl5 = JunctionPSI5VCFDataloader(
                self.combined_splicemap5.to_frame(index=False), fasta_file, vcf_file, encode=False)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl5, event_type='psi5'))

        if self.combined_splicemap3 is not None:
            self.dl3 = JunctionPSI3VCFDataloader(
                self.combined_splicemap3.to_frame(index=False), fasta_file, vcf_file, encode=False)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl3, event_type='psi3'))

    def _iter_dl(self, dl, event_type):
        for row in dl:
            junction_id = row['metadata']['exon']['junction']
            chrom, start, end, strand = junction_str_to_tuple(junction_id)
            row['metadata']['junction'] = dict()
            row['metadata']['junction']['junction'] = junction_id
            row['metadata']['junction']['event_type'] = event_type
            row['metadata']['junction']['Chromosome'] = chrom
            row['metadata']['junction']['Start'] = start
            row['metadata']['junction']['End'] = end
            row['metadata']['junction']['Strand'] = strand
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
