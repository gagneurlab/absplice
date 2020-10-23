from tqdm import tqdm
import pandas as pd
from splicing_outlier_prediction.utils import get_abs_max_rows


class SplicingOutlierResult:

    def __init__(self, df, df_spliceAI=None):
        self.df = df
        self.df_spliceAI = df_spliceAI
        self._junction = None
        self._splice_site = None
        self._gene = None

    @classmethod
    def read_csv(cls, path, **kwargs):
        return cls(pd.read_csv(path, **kwargs))

    def add_spliceAI(self, df):
        """
        Includes spliceAI predictions into results.

        Args:
          df: path to csv or dataframe of spliceAI predictions
        """
        self._gene = None
        if type(df) == str:
            df = pd.read_csv(df)
        self.df_spliceAI = df

    @staticmethod
    def _explode_samples(df):
        df = df.copy()
        df['samples'] = df['samples'].str.split(';')
        return df.rename(columns={'samples': 'sample'}).explode('sample')

    @property
    def psi5(self):
        return SplicingOutlierResult(self.df[self.df['event_type'] == 'psi5'])

    @property
    def psi3(self):
        return SplicingOutlierResult(self.df[self.df['event_type'] == 'psi3'])

    @property
    def junction(self):
        if self._junction is None:
            df = self.df
            if 'samples' in self.df:
                df = self._explode_samples(df)
                index = ['junction', 'sample', 'event_type']
            else:
                index = ['junction', 'event_type']
            self._junction = get_abs_max_rows(
                df.set_index('junction'), index, 'delta_psi') \
                .reset_index('event_type')
        return self._junction

    @property
    def splice_site(self):
        if self._splice_site is None:
            index = ['splice_site', 'event_type']
            if 'samples' in self.df:
                index.append('sample')
            self._splice_site = get_abs_max_rows(
                self.junction, index, 'delta_psi') \
                .reset_index('event_type')
        return self._splice_site

    @property
    def gene(self):
        if self._gene is None:
            index = ['gene_name']
            if 'samples' in self.df:
                index.append('sample')
            self._gene = get_abs_max_rows(
                self.junction, index, 'delta_psi')

            if self.df_spliceAI is not None:
                df_spliceAI = self.df_spliceAI
                if 'samples' in self.df:
                    df_spliceAI = self._explode_samples(df_spliceAI)

                    df_spliceAI = get_abs_max_rows(
                        df_spliceAI, index, 'delta_score')

                self._gene = self._gene.join(df_spliceAI, how='outer',
                                             rsuffix='_spliceAI')

        return self._gene

    def infer_cat(self, cat_inference, progress=False):
        if 'samples' not in self.df.columns:
            raise ValueError('"samples" column is missing.')

        rows = self.junction.iterrows()
        if progress:
            rows = tqdm(rows, total=self.junction.shape[0])

        df = pd.DataFrame([
            cat_inference.infer(junction, sample, row['event_type'])
            for (junction, sample), row in rows
            if cat_inference.contains(junction, sample, row['event_type'])
        ]).set_index(['junction', 'sample'])
        self._junction = self.junction.join(df)
        self._splice_site = None
        self._gene = None

    @staticmethod
    def _add_maf(df, population, default=-1):
        df['maf'] = df['variant'].map(
            lambda x: population.get(x, default))
        return df

    def add_maf(self, population, default=-1):
        self.df = self._add_maf(self.df, population, default)

    def filter_maf(self, max_num_sample=2, population=None,
                   maf_cutoff=0.001, default=-1):
        if max_num_sample:
            df = self.df[self.df['samples'].str.split(';')
                         .map(len) <= max_num_sample]

        if population:
            df = self._add_maf(df, population, default)
            df = df[df['maf'] <= maf_cutoff]

        return SplicingOutlierResult(df)
