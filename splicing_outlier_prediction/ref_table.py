import pandas as pd


class SplicingRefTable:

    def __init__(self, df_ref):
        self.df = df_ref.set_index('junctions')
        self.method = self._infer_method(self.df)

    @classmethod
    def read_csv(cls, path, **kwargs):
        return cls(pd.read_csv(path, **kwargs))

    @staticmethod
    def _infer_method(df):
        if 'k' in df.columns and 'n' in df.columns:
            return 'kn'
        elif 'alpha' in df.columns and 'beta' in df.columns:
            return 'bb'

    @staticmethod
    def valid():
        raise NotImplementedError()

    @staticmethod
    def download(tissue_name):
        raise NotImplementedError()

    @staticmethod
    def fetch(tissue_name):
        raise NotImplementedError()
