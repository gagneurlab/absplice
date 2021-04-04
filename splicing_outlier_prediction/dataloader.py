from splicing_outlier_prediction import SplicingRefTable
import pandas as pd


class RefTableMixin:

    def __init__(self, ref_tables5=list(), ref_tables3=list(), 
                intron_annotation5=None, 
                intron_annotation3=None, 
                **kwargs):
        if not ref_tables5 and not ref_tables3:
            raise ValueError(
                '`ref_tables5` and `ref_tables3` cannot be both empty')

        self.intron_annotation5 = None
        self.intron_annotation3 = None
        if ref_tables5:
            srf5 = self._read_ref_tables(ref_tables5)
            self.ref_tables5 = srf5.ref_tables
            self.intron_annotation5 = srf5
            if len(ref_tables5) > 1:
                self.intron_annotation5_path = intron_annotation5
                srf5.save_concat_ref_tables(save_path=intron_annotation5)
            else:
                self.intron_annotation5_path = ref_tables5[0]
            
        if ref_tables3:
            srf3 = self._read_ref_tables(ref_tables3)
            self.ref_tables3 = srf3.ref_tables
            self.intron_annotation3 = srf3
            if len(ref_tables3) > 1:
                self.intron_annotation3_path = intron_annotation3
                srf3.save_concat_ref_tables(save_path=intron_annotation3)
            else:
                self.intron_annotation3_path = ref_tables3[0]
            
    @staticmethod
    def _read_ref_tables(path):
        if type(path) is list:
            return SplicingRefTable.read_csv(path)
        elif type(path) == SplicingRefTable:
            return path
        else:
            print(type(path))
            raise ValueError(
                'ref_tables should be path to ref_tables file'
                ' or `SplicingRefTable` object')