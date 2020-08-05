import pandas as pd


class SplicingRefTable:

    def __init__(self, df_ref):
        self.ref_table = df_ref
        self.method = self._infer_method(self.ref_table)

    @staticmethod
    def read_csv(path, **kwargs):
        return pd.read_csv(path, **kwargs)
    
    @staticmethod
    def valid():
        raise NotImplementedError()
    
    @staticmethod
    def download(tissue_name):
        raise NotImplementedError()

    @staticmethod
    def fetch(tissue_name):
        raise NotImplementedError()
        
    @staticmethod
    def _infer_method(ref_table):
        return 'kn'  # 'bb'

    def delta_psi(delta_logit_psi):
        raise NotImplementedError()

    def z_score(delta_logit_psi):
        raise NotImplementedError()

    def pvalue(delta_logit_psi):
        raise NotImplementedError()
