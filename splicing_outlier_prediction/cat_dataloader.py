from splicemap import SpliceCountTable as CountTable
from splicing_outlier_prediction.dataloader import SpliceMapMixin
from splicing_outlier_prediction.utils import delta_logit_PSI_to_delta_PSI, logit
import pandas as pd
import re
from typing import List
from splicemap.splice_map import SpliceMap

class CatInference(SpliceMapMixin):

    def __init__(
        self, count_cat, splicemap_cat5=None, splicemap_cat3=None, 
        splicemap5=None, splicemap3=None, sample_mapping=None, name=None
        ):
        SpliceMapMixin.__init__(self, splicemap5, splicemap3)
        self.ct = self._read_cat_count_table(count_cat, name)
        self.contains_chr = self._contains_chr()
        if sample_mapping:
            self._update_samples(sample_mapping)
        self.samples = set(self.ct.samples)
        self.common_junctions5 = list()
        self.common_junctions3 = list()
        self.tissues5 = list()
        self.tissues3 = list()
        self.splicemap_cat5_provided = False
        self.splicemap_cat3_provided = False

        if self.combined_splicemap5 is not None:
            self.tissues5 = [sm.name for sm in self.splicemaps5]
            self.common_junctions5 = [set(sm5.df['junctions'])\
                .intersection(self.ct.junctions) for sm5 in self.splicemaps5]
            self.ct_cat5 = self.ct.filter_event5(list(set.union(*self.common_junctions5)))
            if splicemap_cat5 is not None:
                self.splicemap_cat5_provided = True
                self.splicemap5_cat = SpliceMapMixin._read_splicemap(splicemap_cat5)[0]
            self.ref_psi5_cat = self.ct_cat5.ref_psi5(annotation=False)
        if self.combined_splicemap3 is not None:
            self.tissues3 = [sm.name for sm in self.splicemaps3]
            self.common_junctions3 = [set(sm3.df['junctions'])\
                .intersection(self.ct.junctions) for sm3 in self.splicemaps3]
            self.ct_cat3 = self.ct.filter_event3(list(set.union(*self.common_junctions3)))
            if splicemap_cat3 is not None:
                self.splicemap_cat3_provided = True
                self.splicemap3_cat = SpliceMapMixin._read_splicemap(splicemap_cat3)[0]
            self.ref_psi3_cat = self.ct_cat3.ref_psi3(annotation=False)

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
            
    def _contains_chr(self):
        return 'chr' in self.ct.junctions[0]

    def _update_samples(self, sample_mapping):
        self.ct.update_samples(sample_mapping)

    def contains(self, junction_id, gene_id, tissue, sample, event_type):
        if sample not in self.samples:
            return False
        else:
            if event_type == 'psi5':
                df = self.splicemaps5[self.tissues5.index(tissue)].df
                if tissue not in self.tissues5:
                    raise ValueError('"CatInference and MMSplice+SpliceMap should use the same SpliceMaps"')
                else:
                    return (df[(df['junctions'] == junction_id) & (df['gene_id'] == gene_id)].shape[0] > 0) \
                        and (junction_id in self.common_junctions5[self.tissues5.index(tissue)])
            elif event_type == 'psi3':
                df = self.splicemaps3[self.tissues3.index(tissue)].df
                if tissue not in self.tissues5:
                    raise ValueError('"CatInference and MMSplice+SpliceMap should use the same SpliceMaps"')
                else:
                    return (df[(df['junctions'] == junction_id) & (df['gene_id'] == gene_id)].shape[0] > 0) \
                        and (junction_id in self.common_junctions3[self.tissues3.index(tissue)])
            else:
                raise ValueError('"event_type" should be "psi5" or "psi3"')

    def infer(self, junction_id, gene_id, tissue, sample, event_type, clip_threshold=0.01):
        if not self.contains(junction_id, gene_id, tissue, sample, event_type):
            return {
                'junction': junction_id,
                'gene_id': gene_id,
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

        splicemap_cat_provided = False

        if event_type == 'psi5':
            ct_cat = self.ct_cat5
            psi_cat_df = ct_cat.psi5
            ref_psi_cat_df = self.ref_psi5_cat.df
            if self.splicemap_cat5_provided:
                splicemap_cat_provided = True
                splicemap_cat_df = self.splicemap5_cat.df.set_index(['junctions', 'gene_id'])
            tissue_cat = ct_cat.name
            splicemap_target_df = self.splicemaps5[self.tissues5.index(tissue)]\
                .df.set_index(['junctions', 'gene_id'])
        elif event_type == 'psi3':
            ct_cat = self.ct_cat3
            psi_cat_df = ct_cat.psi3
            ref_psi_cat_df = self.ref_psi3_cat.df
            if self.splicemap_cat3_provided:
                splicemap_cat_provided = True
                splicemap_cat_df = self.splicemap3_cat.df.set_index(['junctions', 'gene_id'])
            tissue_cat = ct_cat.name
            splicemap_target_df = self.splicemaps3[self.tissues3.index(tissue)]\
                .df.set_index(['junctions', 'gene_id'])
        else:
            raise ValueError('Site should be "psi5" or "psi3"')

        ref_psi_target = splicemap_target_df.loc[(junction_id, gene_id)]['ref_psi'] #TODO: splicemap has junction, gene as unique index now!!!
        psi_cat = psi_cat_df.loc[junction_id, sample]
        ref_psi_cat = ref_psi_cat_df.loc[junction_id]['ref_psi']
        median_n_cat = ref_psi_cat_df.loc[junction_id]['median_n']
        # If junction is in SpliceMap of CAT and there is more statistical power, use SpliceMap
        if splicemap_cat_provided:
            if junction_id in splicemap_cat_df.index:
                if splicemap_cat_df.loc[(junction_id, gene_id)]['n'] > ref_psi_cat_df.loc[junction_id]['n']:
                    ref_psi_cat = splicemap_cat_df.loc[(junction_id, gene_id)]['ref_psi']
                    median_n_cat = splicemap_cat_df.loc[(junction_id, gene_id)]['median_n']
            
        delta_logit_psi_cat = logit(psi_cat, clip_threshold) - \
            logit(ref_psi_cat, clip_threshold)
              
        result_infer = {
            'junction': junction_id,
            'gene_id': gene_id,
            'sample': sample,
            'tissue': tissue,
            'count_cat': ct_cat.df.loc[junction_id, sample],
            'psi_cat': psi_cat,
            'ref_psi_cat': ref_psi_cat,
            'k_cat': ref_psi_cat_df.loc[junction_id]['k'],
            'n_cat': ref_psi_cat_df.loc[junction_id]['n'],
            'median_n_cat': median_n_cat,
            'delta_logit_psi_cat': delta_logit_psi_cat,
            'delta_psi_cat': delta_logit_PSI_to_delta_PSI(
                delta_logit_psi_cat, ref_psi_target, clip_threshold=clip_threshold),
            'tissue_cat': tissue_cat
        }
        assert pd.DataFrame(result_infer, index=[0]).shape[0] == 1
        result_infer = pd.DataFrame(result_infer, index=[0]).astype({
            'junction': pd.StringDtype(),
            'gene_id': pd.StringDtype(),
            'sample': pd.StringDtype(),
            'tissue': pd.StringDtype(),
            'count_cat': 'Int64',
            'psi_cat': 'float64',
            'ref_psi_cat': 'float64',
            'k_cat': 'Int64',
            'n_cat': 'Int64',
            'median_n_cat': 'float64',
            'delta_logit_psi_cat': 'float64',
            'delta_psi_cat': 'float64',
            'tissue_cat': pd.StringDtype(),
        }).to_dict('records')[0]
        
        return result_infer


    # def infer_all(self, event_type):

    #     if event_type == 'psi5':
    #         common_junctions = self.common_junctions5
    #         tissues = self.tissues5
    #     elif event_type == 'psi3':
    #         common_junctions = self.common_junctions3
    #         tissues = self.tissues3

    #     infer_rows = list()
    #     for tissue in tissues:
    #         for sample in self.samples:
    #             for junction in common_junctions[tissues.index(tissue)]:
    #                 if self.contains(junction, tissue, sample, event_type):
    #                     infer_rows.append(
    #                         self.infer(junction, tissue, sample, event_type))
    #     df = pd.DataFrame(infer_rows)
    #     df = df.drop_duplicates().set_index(['junction', 'tissue', 'sample'])

    #     return df