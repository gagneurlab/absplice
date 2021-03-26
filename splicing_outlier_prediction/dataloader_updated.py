from splicing_outlier_prediction import SplicingRefTable
import pandas as pd


class RefTableMixin:

    def __init__(self, ref_tables5=list(), ref_tables3=list(), **kwargs):
        if not ref_tables5 and not ref_tables3:
            raise ValueError(
                '`ref_tables5` and `ref_tables3` cannot be both empty')
        self.ref_tables5 = ref_tables5
        self.ref_tables3 = ref_tables3
        self.intron_annotation5 = None
        self.intron_annotation3 = None
        columns = ['junctions', 'Chromosome', 'Start', 'End', 'Strand']
        if ref_tables5:
            self.intron_annotation5 = pd.concat(
                [self._read_ref_tables(ref_table)[columns] for ref_table in ref_tables5]
            ).drop_duplicates(subset='junctions')

        if ref_tables3:
            self.intron_annotation3 = pd.concat(
                [self._read_ref_tables(ref_table)[columns] for ref_table in ref_tables3]
            ).drop_duplicates(subset='junctions')

    @staticmethod
    def _read_ref_tables(path):
        if type(path) == str:
            return SplicingRefTable.read_csv(path)
        elif type(path) == SplicingRefTable:
            return path
        else:
            print(type(path))
            raise ValueError(
                'ref_tables should be path to ref_tables file'
                ' or `SplicingRefTable` object')
