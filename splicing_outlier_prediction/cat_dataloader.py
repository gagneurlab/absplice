from splicemap import SpliceCountTable as CountTable
from splicing_outlier_prediction.dataloader import SpliceMapMixin
from splicing_outlier_prediction.utils import delta_logit_PSI_to_delta_PSI, logit
import pandas as pd
import re
from typing import List

class CatInference(SpliceMapMixin):

    def __init__(self, count_cat, splicemap5=None, splicemap3=None, sample_mapping=None, name=None):
        SpliceMapMixin.__init__(self, splicemap5, splicemap3)
        self.ct = self._read_cat_count_table(count_cat, name)
        if sample_mapping:
            self._update_samples(sample_mapping)
        self.samples = set(self.ct.samples)
        self.common_junctions5 = list()
        self.common_junctions3 = list()
        self.tissues5 = list()
        self.tissues3 = list()

        if self.combined_splicemap5 is not None:
            self.common_junctions5 = [set(sm5.df['junctions'])\
                .intersection(self.ct.junctions) for sm5 in self.splicemaps5]
            self.ct_cat5 = self.ct.filter_event5(list(set.union(*self.common_junctions5)))
            self.ref_psi5_cat = self.ct_cat5.ref_psi5(annotation=False)
            self.tissues5 = [sm.name for sm in self.splicemaps5]
        if self.combined_splicemap3 is not None:
            self.common_junctions3 = [set(sm3.df['junctions'])\
                .intersection(self.ct.junctions) for sm3 in self.splicemaps3]
            self.ct_cat3 = self.ct.filter_event3(list(set.union(*self.common_junctions3)))
            self.ref_psi3_cat = self.ct_cat3.ref_psi3(annotation=False)
            self.tissues3 = [sm.name for sm in self.splicemaps3]

    @staticmethod
    def _read_cat_count_table(path, name):
        if type(path) is str:
            return CountTable.read_csv(path, name)
        elif type(path) is CountTable:
            return [path]
        else:
            print(type(path))
            raise ValueError(
                '`count_cat` argument should'
                ' be path to cat SpliceCountTable files'
                ' or `SpliceCountTable` object')

    def _update_samples(self, sample_mapping):
        self.ct.update_samples(sample_mapping)

    def contains(self, junction_id, sample, tissue, event_type):
        if sample not in self.samples:
            return False
        else:
            if event_type == 'psi5':
                return junction_id in self.common_junctions5[self.tissues5.index(tissue)]
            elif event_type == 'psi3':
                return junction_id in self.common_junctions3[self.tissues3.index(tissue)]
            else:
                raise ValueError('"event_type" should be "psi5" or "psi3"')

    def infer(self, junction_id, sample, tissue, event_type, clip_threshold=0.01):
        if not self.contains(junction_id, sample, tissue, event_type):
            return {
                'junction': junction_id,
                'sample': sample,
                'tissue': tissue,
                'count_cat': None,
                'psi_cat': None,
                'ref_psi_cat': None,
                'k_cat': None,
                'n_cat': None,
                'delta_logit_psi_cat': None,
                'delta_psi_cat': None,
                'tissue_cat': self.ct.name
            }

        if event_type == 'psi5':
            ct_cat = self.ct_cat5
            psi_cat_df = ct_cat.psi5
            ref_psi_cat_df = self.ref_psi5_cat.df
            tissue_cat = ct_cat.name
            splicemap_target_df = self.splicemaps5[self.tissues5.index(tissue)]\
                .df.set_index('junctions')
        elif event_type == 'psi3':
            ct_cat = self.ct_cat3
            psi_cat_df = ct_cat.psi3
            ref_psi_cat_df = self.ref_psi3_cat.df
            tissue_cat = ct_cat.name
            splicemap_target_df = self.splicemaps3[self.tissues3.index(tissue)]\
                .df.set_index('junctions')
        else:
            raise ValueError('Site should be "psi5" or "psi3"')

        psi_cat = psi_cat_df.loc[junction_id, sample]
        ref_psi_cat = ref_psi_cat_df.loc[junction_id]['ref_psi']
        ref_psi_target = splicemap_target_df.loc[junction_id]['ref_psi']
        delta_logit_psi_cat = logit(psi_cat, clip_threshold) - \
            logit(ref_psi_cat, clip_threshold)

        return {
            'junction': junction_id,
            'sample': sample,
            'tissue': tissue,
            'count_cat': ct_cat.df.loc[junction_id, sample],
            'psi_cat': psi_cat,
            'ref_psi_cat': ref_psi_cat,
            'k_cat': ref_psi_cat_df.loc[junction_id]['k'],
            'n_cat': ref_psi_cat_df.loc[junction_id]['n'],
            'median_n_cat': ref_psi_cat_df.loc[junction_id]['median_n'],
            'delta_logit_psi_cat': delta_logit_psi_cat,
            'delta_psi_cat': delta_logit_PSI_to_delta_PSI(
                delta_logit_psi_cat, ref_psi_target, clip_threshold=clip_threshold),
            'tissue_cat': tissue_cat
        }
