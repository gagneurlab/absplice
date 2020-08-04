import pandas as pd


class SplicingRefTable:

    def __init__(self, ref_table):
        self.ref_table = pd.read_csv(ref_table)
        self.method = self._infer_method(ref_table)

    @staticmethod
    def _infer_method(ref_table):
        return 'kn'  # 'bb'

    def delta_psi(delta_logit_psi):
        raise NotImplementedError()

    def z_score(delta_logit_psi):
        raise NotImplementedError()

    def pvalue(delta_logit_psi):
        raise NotImplementedError()
