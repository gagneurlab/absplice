import pandas as pd


class SplicingRefTable:

    def __init__(self, df_ref_tables):
        self.ref_tables = df_ref_tables
        self.df = self._concat_ref_tables(self.ref_tables)
        self.method = self._infer_method(self.df)
        
    @classmethod
    def read_csv(cls, path_list, **kwargs):
        l = list()
        [l.append(pd.read_csv(path, **kwargs)) for path in path_list]
        return cls(l)

    @staticmethod
    def _concat_ref_tables(ref_tables):
        # columns = ['junctions', 'Chromosome', 'Start', 'End', 'Strand']
        df = pd.concat(
            [ref_table for ref_table in ref_tables]
            ).drop_duplicates(subset='junctions').set_index('junctions')
        return df

    def save_concat_ref_tables(self, save_path):
        columns = ['junctions', 'Chromosome', 'Start', 'End', 'Strand']
        self.df.reset_index()[columns].to_csv(save_path, index=False)

    @staticmethod
    def _infer_method(df):
        if 'k' in df.columns and 'n' in df.columns:
            return 'kn'
        elif 'alpha' in df.columns and 'beta' in df.columns:
            return 'bb'

    @property
    def junctions(self):
        return self.df.index

    @staticmethod
    def valid():
        raise NotImplementedError()

    @staticmethod
    def download(tissue_name):
        raise NotImplementedError()

    @staticmethod
    def fetch(tissue_name):
        raise NotImplementedError()
