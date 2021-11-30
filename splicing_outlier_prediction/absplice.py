from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
from splicing_outlier_prediction.utils import get_abs_max_rows
from splicing_outlier_prediction.cat_dataloader import CatInference


class AbSplice:

    def __init__(self, df_mmsplice, df_spliceAI)
        self.df_mmsplice = df_mmsplice
        self.df_spliceAI = df_spliceAI
        self.df_absplice = self._join(df_mmsplice, df_spliceAI)
        self._variant = None
        self._gene = None


    @classmethod
    def read_csv(cls, path, **kwargs):
        return cls(pd.read_csv(path, **kwargs))

    
    @staticmethod
    def _join(df_mmsplice, df_spliceAI):
        return df_mmsplice.join(df_spliceAI, on=['variant', 'gene_id'], how='outer')


    @staticmethod
    def _explode(df, col='samples', new_name='sample'):
        if new_name == None:
            new_name = col
        df = df.copy()
        df[col] = df[col].str.split(';').map(
            set, na_action='ignore').map(list, na_action='ignore')
        return df.rename(columns={col: new_name}).explode(new_name)



    @property
    def variant(self, df):
        """
        max aggregate scores on variant level (if provided, for each sample and tissue)
        stores scores in self._variant
        """
        if self._variant is None:
            df = self.df
            if 'samples' in self.df:
                df = self._explode(df, col='samples', new_name='sample')
                index = ['variant', 'sample']
            else:
                index = ['variant']
            if 'tissue' in self.df:
                index.append('tissue')
            self._variant = get_abs_max_rows(
                df.set_index('variant'), index, 'delta_psi') \
                .reset_index('event_type')
        return self._variant