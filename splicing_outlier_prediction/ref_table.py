import pandas as pd
import re


class SplicingRefTable:

    def __init__(self, df_ref_tables):
        self.ref_tables = df_ref_tables
        self.df = self._combine_ref_tables(self.ref_tables)
        self.df_all = self._concat_ref_tables(self.ref_tables)
        self.method = self._infer_method(self.ref_tables)
        
    @classmethod
    def read_csv(cls, path_list, regex_pattern=None, **kwargs):
        l = list()
        [l.append(pd.read_csv(path, **kwargs)) for path in path_list]
        if regex_pattern:
            for (i, path) in enumerate(path_list):
                if 'tissue' not in l[i].columns:
                    if isinstance(regex_pattern, str):
                        tissue = re.search(regex_pattern, path).group(1)
                    elif isinstance(regex_pattern, list):
                        tissue = regex_pattern[i]
                    # tissue = re.search(regex_pattern, path).group(1)
                    l[i]['tissue'] = tissue
        return cls(l)

    @staticmethod
    def _combine_ref_tables(ref_tables):
        columns = ['junctions', 'Chromosome', 'Start', 'End', 'Strand', ]
        df = pd.concat(
            [ref_table[columns] for ref_table in ref_tables]
            ).drop_duplicates(subset='junctions').set_index('junctions')
        return df
    
    @staticmethod
    def _concat_ref_tables(ref_tables):
        df = pd.concat(
            [ref_table for ref_table in ref_tables]
            ).set_index('junctions')
        return df

    def save_combined_ref_tables(self, save_path):
        self.df.reset_index().to_csv(save_path, index=False)

    @staticmethod
    def _infer_method(ref_tables):
        method = list()
        for df in ref_tables:
            if 'k' in df.columns and 'n' in df.columns:
                method.append('kn')
            elif 'alpha' in df.columns and 'beta' in df.columns:
                method.append('bb')
        return method

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
