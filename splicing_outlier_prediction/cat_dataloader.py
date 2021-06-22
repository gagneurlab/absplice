 

from splicemap import SpliceCountTable as CountTable
from splicing_outlier_prediction.dataloader import SpliceMapMixin
from splicing_outlier_prediction.utils import delta_logit_PSI_to_delta_PSI, logit
import pandas as pd
import re
from typing import List


class CatInferenceMixin(SpliceMapMixin):

    def __init__(self, count_cat, splicemap5=None, splicemap3=None):
        SpliceMapMixin.__init__(self, splicemap5, splicemap3)
        self.ct = self._read_cat_count_table(count_cat)
        self.samples = self._samples(self.ct)
        self.common_junctions5 = list()
        self.common_junctions3 = list()
 
        if self.combined_splicemap5 is not None:
            self.common_junctions5 = self._common_junctions5(self.ct, self.combined_splicemap5)
            self.ct_cat5 = self._ct_cat5(self.ct, self.common_junctions5)
            self.ref_psi5_cat = self._ref_psi5_cat(self.ct_cat5)
        if self.combined_splicemap3 is not None:
            self.common_junctions3 = self._common_junctions3(self.ct, self.combined_splicemap3)
            self.ct_cat3 = self._ct_cat3(self.ct, self.common_junctions3)
            self.ref_psi3_cat = self._ref_psi3_cat(self.ct_cat3)

    @staticmethod
    def _read_cat_count_table(path):
        if type(path) is str:
            return [CountTable.read_csv(path)]
        elif type(path) is CountTable:
            return [path]
        elif type(path) is list:
            return [CatInferenceMixin._read_cat_count_table(i)[0] for i in path]
        else:
            print(type(path))
            raise ValueError(
                '`count_cat` argument should'
                ' be list of path to cat SpliceCountTable files'
                ' or `SpliceCountTable` object')

    @staticmethod
    def _samples(cat_count_tables: List[CountTable]):
        return [set(ct.samples) for ct in cat_count_tables]

    @staticmethod
    def _common_junctions5(cat_count_tables: List[CountTable], combined_splicemap5):
        return [set(combined_splicemap5.index).intersection(ct.junctions) for ct in cat_count_tables]

    @staticmethod
    def _common_junctions3(cat_count_tables: List[CountTable], combined_splicemap3):
        return [set(combined_splicemap3.index).intersection(ct.junctions) for ct in cat_count_tables]

    @staticmethod
    def _ct_cat5(cat_count_tables: List[CountTable], common_junctions5):
        return [ct.filter_event5(common_j5) for (ct, common_j5) in zip(cat_count_tables, common_junctions5)]

    @staticmethod
    def _ct_cat3(cat_count_tables: List[CountTable], common_junctions3):
        return [ct.filter_event3(common_j3) for (ct, common_j3) in zip(cat_count_tables, common_junctions3)]

    @staticmethod
    def _ref_psi5_cat(ct_cat5):
        return [ct.ref_psi5(annotation=False) for ct in ct_cat5]

    @staticmethod
    def _ref_psi3_cat(ct_cat3):
        return [ct.ref_psi3(annotation=False) for ct in ct_cat3]

        
class CatInference(CatInferenceMixin):

    def __init__(self, count_cat, splicemap5=None, splicemap3=None, sample_mapping=None):
        CatInferenceMixin.__init__(self, count_cat, splicemap5, splicemap3)

        if sample_mapping:
            self._update_samples(sample_mapping)
            self.samples = self._samples(self.ct)

    def _update_samples(self, sample_mapping):
        [self.ct[i].update_samples(sample_mapping) for i in range(len(self.ct))]

    def contains(self, junction_id, sample, event_type):
        if sample not in set.union(*self.samples):
            return False
        else:
            if event_type == 'psi5':
                return junction_id in set.union(*self.common_junctions5)
            elif event_type == 'psi3':
                return junction_id in set.union(*self.common_junctions3)
            else:
                raise ValueError('"event_type" should be "psi5" or "psi3"')

    def _which_cat_contains(self, junction_id, sample, event_type):
        if event_type == 'psi5':
            return [junction_id in j and sample in s for (j, s) in zip(self.common_junctions5,self.samples)]
        if event_type == 'psi3':
            return [junction_id in j and sample in s for (j, s) in zip(self.common_junctions3,self.samples)]

    def _which_target_contains(self, junction_id, event_type):
        if event_type == 'psi5':
            return [junction_id in set(sm.df['junctions']) for sm in self.splicemaps5]
        if event_type == 'psi3':
            return [junction_id in set(sm.df['junctions']) for sm in self.splicemaps3]

    def _infer_single(self, i, j, junction_id, sample, event_type, clip_threshold):
        if event_type == 'psi5':
            ct_cat = self.ct_cat5[j]
            psi_cat_df = ct_cat.psi5
            ref_psi_cat_df = self.ref_psi5_cat[j].df
            tissue_cat = self.ct_cat5[j].name
            splicemap_target_df = self.splicemaps5[i].df.set_index('junctions')
            tissue_target = self.splicemaps5[i].name
        elif event_type == 'psi3':
            ct_cat = self.ct_cat3[j]
            psi_cat_df = ct_cat.psi3
            ref_psi_cat_df = self.ref_psi3_cat[j].df
            tissue_cat = self.ct_cat3[j].name
            splicemap_target_df = self.splicemaps3[i].df.set_index('junctions')
            tissue_target = self.splicemaps3[i].name
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
            'tissue': tissue_target,
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

    def infer(self, junction_id, sample, event_type, clip_threshold=0.01):
        cats_with_junction = self._which_cat_contains(junction_id, sample, event_type)
        targets_with_junction = self._which_target_contains(junction_id, event_type)

        cat_infer_results = list()
        for i in range(len(targets_with_junction)):
            for j in range(len(cats_with_junction)):
                if targets_with_junction[i] and cats_with_junction[j]:
                    cat_infer_single = self._infer_single(i, j, junction_id, sample, event_type, clip_threshold)
                    cat_infer_results.append(cat_infer_single)

        return cat_infer_results



    
