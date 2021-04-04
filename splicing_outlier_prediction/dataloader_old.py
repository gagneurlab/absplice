from splicing_outlier_prediction import SplicingRefTable


class RefTableMixin:

    def __init__(self, ref_table5=None, ref_table3=None, **kwargs):
        if ref_table5 is None and ref_table3 is None:
            raise ValueError(
                '`ref_table5` and `ref_table3` cannot be both None')
        self.ref_table5 = None
        self.ref_table3 = None
        if ref_table5:
            self.ref_table5 = self._read_ref_table(ref_table5)
        if ref_table3:
            self.ref_table3 = self._read_ref_table(ref_table3)

    @staticmethod
    def _read_ref_table(path):
        if type(path) == str:
            return SplicingRefTable.read_csv(path)
        elif type(path) == SplicingRefTable:
            return path
        else:
            raise ValueError(
                'ref_table should be path to ref_table file'
                ' or `SplicingRefTable` object')
