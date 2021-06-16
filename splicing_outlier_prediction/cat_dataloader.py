from splicemap import SpliceCountTable as CountTable
from splicing_outlier_prediction.dataloader import SpliceMapMixin
from splicing_outlier_prediction.utils import delta_logit_PSI_to_delta_PSI, logit
import pandas as pd
import re


class CatInference(SpliceMapMixin):

    def __init__(self, count_cat:  CountTable, splicemap5=None, splicemap3=None):
        SpliceMapMixin.__init__(self, splicemap5, splicemap3)
        self.ct = list()
        self.samples = list()
        self.common_junctions5 = list()
        self.common_junctions3 = list()
        self.ct_cat5 = list()
        self.ct_cat3 = list()
        self.ref_psi5_cat = list()
        self.ref_psi3_cat = list()
        self.tissue_cat5 = list()
        self.tissue_cat3 = list()

        for i in range(len(count_cat)):
            self.ct.append(CountTable.read_csv(count_cat[i]))
            if sample_mapping is not None:
                self.ct[i].update_samples(sample_mapping)
            self.samples.append(set(self.ct[i].samples))
            if self.ref_tables5 is not None:
                self.common_junctions5.append(set(self.combined_ref_tables5.df_all.index)
                                              .intersection(self.ct[i].junctions))
                self.ct_cat5.append(self.ct[i].filter_event5(
                    self.common_junctions5[i]))
                self.ref_psi5_cat.append(
                    self.ct_cat5[i].ref_psi5(annotation=False))
                self.tissue_cat5.append(self._annotate_tissue_cat(
                    count_cat, i, regex_pattern_cat))
            if self.ref_tables3 is not None:
                self.common_junctions3.append(set(self.combined_ref_tables3.df_all.index)
                                              .intersection(self.ct[i].junctions))
                self.ct_cat3.append(self.ct[i].filter_event3(
                    self.common_junctions3[i]))
                self.ref_psi3_cat.append(
                    self.ct_cat3[i].ref_psi3(annotation=False))
                self.tissue_cat3.append(self._annotate_tissue_cat(
                    count_cat, i, regex_pattern_cat))

    def _annotate_tissue_cat(self, count_cat, i, regex_pattern_cat):
        tissue_cat = None
        if regex_pattern_cat:
            if isinstance(regex_pattern_cat, str):
                tissue_cat = re.search(
                    regex_pattern_cat, count_cat[i]).group(1)
            elif isinstance(regex_pattern_cat, list):
                tissue_cat = regex_pattern_cat[i]
            else:
                raise(
                    '"regex_pattern_cat" has to be either a list of tissue names or a regular expression')
        return tissue_cat

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

    def _update_cat_infer_results(self):
        return {
            'junction': self._junction_id,
            'sample': self._sample,
            'tissue': self._tissue,
            'count_cat': self._count_cat,
            'psi_cat': self._psi_cat,
            'ref_psi_cat': self._ref_psi_cat,
            'k_cat': self._k_cat,
            'n_cat': self._n_cat,
            'median_n_cat': self._median_n_cat,
            'delta_logit_psi_cat': self._delta_logit_psi_cat,
            'delta_psi_cat': self._delta_psi_cat,
            'tissue_cat': self._tissue_cat
        }

    def infer(self, junction_id, sample, event_type, clip_threshold=0.01):
        # TODO: if junction in multiple ref tables, creates duplicates...
        self._sample = sample
        self._junction_id = junction_id
        cat_infer_results = list()
        for i in range(len(self.ct)):  # loop through all CAT
            if event_type == 'psi5':
                ref_tables = self.combined_ref_tables5.df_all
                if sample in self.samples[i]:
                    ct_cat = self.ct_cat5[i]
                    psi_cat = ct_cat.psi5
                    ref_psi_cat = self.ref_psi5_cat[i]
                    tissue_cat = self.tissue_cat5[i]
                    common_junctions = self.common_junctions5[i]
            elif event_type == 'psi3':
                ref_tables = self.combined_ref_tables3.df_all
                if sample in self.samples[i]:
                    ct_cat = self.ct_cat3[i]
                    psi_cat = ct_cat.psi3
                    ref_psi_cat = self.ref_psi3_cat[i]
                    tissue_cat = self.tissue_cat3[i]
                    common_junctions = self.common_junctions3[i]
            else:
                raise ValueError('Site should be "psi5" or "psi3"')

            # infer delta psi for current CAT
            if sample in self.samples[i] and junction_id in common_junctions:
                psi = psi_cat.loc[junction_id, sample]
                ref_psi = ref_psi_cat.loc[junction_id]['ref_psi']
                ref_psi_target = ref_tables.loc[junction_id]['ref_psi']
                delta_logit_psi = logit(psi, clip_threshold) - \
                    logit(ref_psi, clip_threshold)
                delta_psi_cat_infer = delta_logit_PSI_to_delta_PSI(
                    delta_logit_psi, ref_psi_target, clip_threshold=clip_threshold)
                self._count_cat = ct_cat.df.loc[junction_id, sample]
                self._psi_cat = psi
                self._ref_psi_cat = ref_psi
                self._k_cat = ref_psi_cat.loc[junction_id]['k']
                self._n_cat = ref_psi_cat.loc[junction_id]['n']
                self._median_n_cat = ref_psi_cat.loc[junction_id]['median_n']
                self._delta_logit_psi_cat = delta_logit_psi
                self._tissue_cat = tissue_cat

                tissue_target = ref_tables.loc[junction_id]['tissue']
                # loop through all non-CAT that contain junction_id
                if type(ref_tables.loc[junction_id]['tissue']) == pd.core.series.Series:
                    for j in range(ref_tables.loc[junction_id]['tissue'].shape[0]):
                        self._tissue = tissue_target.iloc[j]
                        if sample in self.samples[i]:
                            self._delta_psi_cat = delta_psi_cat_infer.iloc[j]
                        cat_infer_results.append(
                            self._update_cat_infer_results())
                else:
                    self._tissue = tissue_target
                    if sample in self.samples[i]:
                        self._delta_psi_cat = delta_psi_cat_infer
                    cat_infer_results.append(self._update_cat_infer_results())

        return cat_infer_results
