from collections import defaultdict
import itertools
from typing import List
import pandas as pd
from tqdm import tqdm
from kipoi.data import SampleIterator
from splicemap.splice_map import SpliceMap

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
            self.combined_splicemap5 = self._combine_junctions(
                self.splicemaps5)
            self.metadata_splicemap5 = self._splicemap_metadata(self.splicemaps5)
        else:
            self.combined_splicemap5 = None

        if splicemap3 is not None:
            self.splicemaps3 = self._read_splicemap(splicemap3)
            self.combined_splicemap3 = self._combine_junctions(
                self.splicemaps3)
            self.metadata_splicemap3 = self._splicemap_metadata(self.splicemaps3)
        else:
            self.combined_splicemap3 = None

    def _splicemap_metadata(self, splicemaps):
        metadata = defaultdict(list)
        
        cols = ['junction', 'gene_id', 'tissue', 'ref_psi', 'median_n', 'gene_name', 'splice_site']
            
        for splicemap in splicemaps:
            df = splicemap.df.copy()
            df = df.rename(columns={'junctions': 'junction'})
            df['tissue'] = splicemap.name 
            itertuples = df[cols].itertuples(index=False)
            
            if self.progress:
                itertuples = tqdm(itertuples)
            
            for row in itertuples:
                metadata[row.junction].append(tuple(row))       
        
        return dict(metadata)
    
    @staticmethod
    def _combine_junctions(splicemaps: List[SpliceMap]):
        columns = ['junctions', 'Chromosome', 'Start', 'End', 'Strand']
        df = pd.concat(
            [s.df[columns] for s in splicemaps]
        ).drop_duplicates(subset='junctions').set_index('junctions')
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
                self.combined_splicemap5, fasta_file, vcf_file, encode=False)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl5, self.combined_splicemap5, event_type='psi5'))

        if self.combined_splicemap3 is not None:
            self.dl3 = JunctionPSI3VCFDataloader(
                self.combined_splicemap3, fasta_file, vcf_file, encode=False)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl3, self.combined_splicemap3, event_type='psi3'))

    def _iter_dl(self, dl, intron_annotations, event_type):
        for row in dl:
            junction_id = row['metadata']['exon']['junction']
            ref_row = intron_annotations.loc[junction_id]
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
