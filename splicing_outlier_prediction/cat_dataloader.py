from kipoi.data import Dataset
from count_table import CountTable
from splicing_outlier_prediction import SplicingRefTable
from splicing_outlier_prediction.dataloader import RefTableMixin


class CatDataloader(RefTableMixin, Dataset):

    def __init__(self, count_cat, ref_table5=None, ref_table3=None, **kwargs):
        super().__init__(ref_table5, ref_table3, **kwargs)

        # TODO: bug in count table because columns[:4] is wrong is index in columns
        self.ct = CountTable.read_csv(count_cat, **kwargs)
        self.samples = self.ct.samples

        self.common_junctions5 = list()
        self.common_junctions3 = list()

        if self.ref_table5:
            self.common_junctions5 = list(
                set(self.ref_table5.junctions).intersection(self.ct.junctions))
            self.count_cat5 = self.ct.filter_event5(self.common_junctions5)
            self.ref_psi5_cat = self.count_cat5.ref_psi5(annotation=False)
        if self.ref_table3:
            self.common_junctions3 = list(
                set(self.ref_table3.junctions).intersection(self.ct.junctions))
            self.count_cat3 = self.ct.filter_event3(self.common_junctions3)
            self.ref_psi3_cat = self.count_cat3.ref_psi3(annotation=False)

    def __len__(self):
        return len(self.common_junctions5) + len(self.common_junctions3)

    def __getitem__(self, idx):
        if idx < len(self.common_junctions5):
            event_type = 'psi5'
            junction_id = self.common_junctions5[idx]
            count_cat = self.count_cat5
            psi_cat = count_cat.psi5
            ref_psi_cat = self.ref_psi5_cat
            ref_table = self.ref_table5
        elif idx < len(self):
            idx = idx - len(self.common_junctions5)
            event_type = 'psi3'
            junction_id = self.common_junctions3[idx]
            count_cat = self.count_cat3
            psi_cat = count_cat.psi3
            ref_psi_cat = self.ref_psi3_cat
            ref_table = self.ref_table3
        else:
            raise IndexError

        return {
            "inputs": {
                'counts': count_cat.df.loc[junction_id, self.samples].values,
                'psi': psi_cat.loc[junction_id, self.samples].values
            },
            "metadata": {
                "junction": junction_id,
                "event_type": event_type,
                'cat_tissue': {
                    'psi_ref': ref_psi_cat.loc[junction_id]['ref_psi'],
                    'k': ref_psi_cat.loc[junction_id]['k'],
                    'n': ref_psi_cat.loc[junction_id]['n'],
                },
                'target_tissue': ref_table.df.loc[junction_id].to_dict()
            }
        }
