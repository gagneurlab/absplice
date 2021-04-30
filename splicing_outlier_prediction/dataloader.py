from splicing_outlier_prediction import SplicingRefTable
import pandas as pd


class RefTableMixin:

    def __init__(self, ref_tables5=list(), ref_tables3=list(), 
                combined_ref_tables5=None, 
                combined_ref_tables3=None, 
                regex_pattern=None,
                save_combined_ref_tables=True,
                **kwargs):
        if not ref_tables5 and not ref_tables3:
            raise ValueError(
                '`ref_tables5` and `ref_tables3` cannot be both empty')

        self.combined_ref_tables5 = None
        self.combined_ref_tables3 = None
        self.ref_tables5 = None
        self.ref_tables3 = None

        if ref_tables5:
            srf5 = self._read_ref_tables(ref_tables5, regex_pattern)
            self.ref_tables5 = srf5.ref_tables
            self.combined_ref_tables5 = srf5
            if len(ref_tables5) > 1:
                self.combined_ref_tables5_path = combined_ref_tables5
                if save_combined_ref_tables:
                    srf5.save_combined_ref_tables(save_path=combined_ref_tables5)
            else:
                self.combined_ref_tables5_path = ref_tables5[0]
            
        if ref_tables3:
            srf3 = self._read_ref_tables(ref_tables3, regex_pattern)
            self.ref_tables3 = srf3.ref_tables
            self.combined_ref_tables3 = srf3
            if len(ref_tables3) > 1:
                self.combined_ref_tables3_path = combined_ref_tables3
                if save_combined_ref_tables:
                    srf3.save_combined_ref_tables(save_path=combined_ref_tables3)
            else:
                self.combined_ref_tables3_path = ref_tables3[0]
            
    @staticmethod
    def _read_ref_tables(path, regex_pattern):
        if type(path) is list:
            return SplicingRefTable.read_csv(path, regex_pattern)
        elif type(path) == SplicingRefTable:
            return path
        else:
            print(type(path))
            raise ValueError(
                'ref_tables should be list of path to ref_table files'
                ' or `SplicingRefTable` object')