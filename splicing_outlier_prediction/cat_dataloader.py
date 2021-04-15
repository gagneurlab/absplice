from count_table import CountTable
from splicing_outlier_prediction.dataloader import RefTableMixin
from splicing_outlier_prediction.utils import delta_logit_PSI_to_delta_PSI, logit


class CatInference(RefTableMixin):

    def __init__(self, count_cat, ref_tables5=None, ref_tables3=None, **kwargs):
        super().__init__(ref_tables5, ref_tables3, **kwargs)
        self.ct = CountTable.read_csv(count_cat)
        self.samples = set(self.ct.samples)
        self.common_junctions5 = list()
        self.common_junctions3 = list()

        if self.combined_ref_tables5 is not None:
            self.common_junctions5 = set(self.combined_ref_tables5.junctions) \
                .intersection(self.ct.junctions)
            self.ct_cat5 = self.ct.filter_event5(self.common_junctions5)
            self.ref_psi5_cat = self.ct_cat5.ref_psi5(annotation=False)
        if self.combined_ref_tables3 is not None:
            self.common_junctions3 = set(self.combined_ref_tables3.junctions) \
                .intersection(self.ct.junctions)
            self.ct_cat3 = self.ct.filter_event3(self.common_junctions3)
            self.ref_psi3_cat = self.ct_cat3.ref_psi3(annotation=False)

    def contains(self, junction_id, sample, event_type):
        if sample not in self.samples:
            return False
        else:
            if event_type == 'psi5':
                return junction_id in self.common_junctions5
            elif event_type == 'psi3':
                return junction_id in self.common_junctions3
            else:
                raise ValueError('"event_type" should be "psi5" or "psi3"')

    def infer(self, junction_id, sample, event_type, clip_threshold=0.01):
        if sample not in self.samples:
            return {
                'junction': junction_id,
                'sample': sample,
                'count_cat': None,
                'psi_cat': None,
                'ref_psi_cat': None,
                'k_cat': None,
                'n_cat': None,
                'delta_logit_psi_cat': None,
                'delta_psi_cat': None
            }

        if event_type == 'psi5':
            ct_cat = self.ct_cat5
            psi_cat = ct_cat.psi5
            ref_psi_cat = self.ref_psi5_cat
            ref_tables = self.combined_ref_tables5.df
        elif event_type == 'psi3':
            ct_cat = self.ct_cat3
            psi_cat = ct_cat.psi3
            ref_psi_cat = self.ref_psi3_cat
            ref_tables = self.combined_ref_tables3.df
        else:
            raise ValueError('Site should be "psi5" or "psi3"')

        psi = psi_cat.loc[junction_id, sample]
        ref_psi = ref_psi_cat.loc[junction_id]['ref_psi']
        ref_psi_target = ref_tables.loc[junction_id]['ref_psi']
        delta_logit_psi = logit(psi, clip_threshold) - \
            logit(ref_psi, clip_threshold)

        return {
            'junction': junction_id,
            'sample': sample,
            'count_cat': ct_cat.df.loc[junction_id, sample],
            'psi_cat': psi,
            'ref_psi_cat': ref_psi,
            'k_cat': ref_psi_cat.loc[junction_id]['k'],
            'n_cat': ref_psi_cat.loc[junction_id]['n'],
            'median_n_cat': ref_psi_cat.loc[junction_id]['median_n'],
            'delta_logit_psi_cat': delta_logit_psi,
            'delta_psi_cat': delta_logit_PSI_to_delta_PSI(
                delta_logit_psi, ref_psi_target, clip_threshold=clip_threshold)
        }
