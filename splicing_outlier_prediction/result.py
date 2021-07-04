from tqdm import tqdm
import pandas as pd
import pickle
from splicing_outlier_prediction.utils import get_abs_max_rows


class SplicingOutlierResult:

    def __init__(self, df, df_spliceAI=None):
        self.df = df
        self.df_spliceAI = df_spliceAI
        self._junction = None
        self._splice_site = None
        self._gene = None
        self._junction_cat_concat = None
        self._splice_site_cat_concat = None
        self._gene_cat_concat = None
        self._junction_cat_features = None
        self._splice_site_cat_features = None
        self._gene_cat_features = None

    @classmethod
    def read_csv(cls, path, **kwargs):
        return cls(pd.read_csv(path, **kwargs))

    @staticmethod
    def _explode_samples(df):
        df = df.copy()
        df['samples'] = df['samples'].str.split(';').map(
            set, na_action='ignore').map(list, na_action='ignore')
        return df.rename(columns={'samples': 'sample'}).explode('sample')

    @property
    def psi5(self):
        return SplicingOutlierResult(self.df[self.df['event_type'] == 'psi5'])

    @property
    def psi3(self):
        return SplicingOutlierResult(self.df[self.df['event_type'] == 'psi3'])

    def filter_samples_with_RNA_seq(self, samples_for_tissue):
        if 'tissue' not in self.df or 'samples' not in self.df:
            raise KeyError('tissue or sample annotation missing')

        if self.df is not None:
            l = list()
            for tissue, samples in samples_for_tissue.items():
                df = self.df[self.df['tissue'] == tissue]
                df = df[~df['samples'].isna()]
                df['samples'] \
                    = df['samples'].apply(lambda x: ';'.join(i for i in x.split(';') if i in samples))
                df_filtered = df[(df['tissue'] == tissue)
                                 & ~(df['samples'] == '')]
                l.append(df_filtered)
            self.df = pd.concat(l)

        if self.df_spliceAI is not None:
            if 'tissue' not in self.df_spliceAI:
                df_spliceAI, index_spliceAI = self._add_tissue_info_to_spliceAI()
            else:
                df_spliceAI = self.df_spliceAI.copy()
            l = list()
            for tissue, samples in samples_for_tissue.items():
                df = df_spliceAI[df_spliceAI['tissue'] == tissue]
                df = df[~df['samples'].isna()]
                df['samples'] \
                    = df['samples'].apply(lambda x: ';'.join(i for i in x.split(';') if i in samples))
                df_filtered = df[(df['tissue'] == tissue)
                                 & ~(df['samples'] == '')]
                l.append(df_filtered)
            self.df_spliceAI = pd.concat(l)

    @property
    def junction(self):
        if self._junction is None:
            df = self.df
            if 'samples' in self.df:
                df = self._explode_samples(df)
                index = ['junction', 'sample', 'event_type']
            else:
                index = ['junction', 'event_type']
            if 'tissue' in self.df:
                index.append('tissue')
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
            if 'tissue' in self.df:
                index.append('tissue')
            if 'tissue_cat' in self.junction:
                index.append('tissue_cat')

            self._splice_site = get_abs_max_rows(
                self.junction, index, 'delta_psi') \
                .reset_index('event_type')
            if 'tissue_cat' in index:
                self._splice_site = self._splice_site.reset_index('tissue_cat')
        return self._splice_site

    @property
    def gene(self):
        if self._gene is None:
            index = ['gene_name']
            if 'samples' in self.df:
                index.append('sample')
            if 'tissue' in self.df and 'tissue' not in index:
                index.append('tissue')
            if 'tissue_cat' in self.junction:
                index.append('tissue_cat')
            self._gene = get_abs_max_rows(
                self.junction, index, 'delta_psi')
            if 'tissue_cat' in index:
                self._gene = self._gene.reset_index('tissue_cat')

            if self.df_spliceAI is not None:
                df_spliceAI, index_spliceAI = self._add_tissue_info_to_spliceAI()
                if 'samples' in self.df_spliceAI:
                    index_spliceAI.append('sample')
                    df_spliceAI = self._explode_samples(df_spliceAI)
                self.df_spliceAI = get_abs_max_rows(
                    df_spliceAI, index_spliceAI, 'delta_score')  # TODO: double check that, before it was only for 'samples' in self.df
                self._gene = self._join_spliceAI(self._gene)

        return self._gene

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

    def _add_tissue_info_to_spliceAI(self):
        df_spliceAI = self.df_spliceAI
        index_spliceAI = ['gene_name']
        if 'tissue' in self.df:
            if 'tissue' in self.df_spliceAI:
                index_spliceAI.append('tissue')
            else:
                l = list()
                for tissue in self.df['tissue'].unique():
                    _df = df_spliceAI.copy()
                    _df['tissue'] = tissue
                    l.append(_df)
                df_spliceAI = pd.concat(l)
                index_spliceAI.append('tissue')

        return df_spliceAI, index_spliceAI

    def _join_spliceAI(self, df):
        index_spliceAI = self.df_spliceAI.index.names
        index = df.index.names
        df = df[df.columns.difference(
            ['index', 'variant_spliceAI', 'delta_score', 'acceptor_gain', 'acceptor_loss', 'donor_gain',
             'donor_loss', 'acceptor_gain_position', 'acceptor_loss_positiin',
             'donor_gain_position', 'donor_loss_position', 'GQ', 'DP_ALT']
        )]
        df = df.reset_index().set_index(index_spliceAI)\
            .join(self.df_spliceAI, how='outer', rsuffix='_spliceAI')\
            .reset_index().set_index(index)
        if 'tissue_cat' in index:
            df = df.reset_index('tissue_cat')
        return df

    def infer_cat(self, cat_inference, progress=False):
        """
        cat_inference: List[CatInference] or CatInference
        """
        if 'samples' not in self.df.columns:
            raise ValueError('"samples" column is missing.')

        if type(cat_inference) == CatInference:
            cat_inference = [cat_inference]

        infer_rows = list()

        for cat in cat_inference:

            rows = self.junction.iterrows()
            if progress:
                rows = tqdm(rows, total=self.junction.shape[0])

            for (junction, sample, tissue), row in rows:
                if cat_inference.contains(junction, sample, row['event_type']):
                    for j in cat_inference.infer(junction, sample, row['event_type']):
                        infer_rows.append(j)

        df = pd.DataFrame(infer_rows)

        # df = pd.DataFrame([
        #     j
        #     for (junction, sample, tissue), row in rows
        #     for j in cat_inference.infer(junction, sample, row['event_type'])
        #     if cat_inference.contains(junction, sample, row['event_type'])
        # ])
        df = df.drop_duplicates().set_index(['junction', 'sample', 'tissue'])
        self._junction = self.junction.join(df)
        self._splice_site = None
        self._gene = None

    @property
    def junction_cat_concat(self):
        if self._junction_cat_concat is None:
            if self._junction is None:
                self.junction
            if 'tissue_cat' not in self._junction:
                raise IndexError(
                    '"tissue_cat not in columns. You have to run infer_cat first"')

            index = ['junction', 'event_type']
            if 'samples' in self.df:
                index.append('sample')
            if 'tissue' in self.df and 'tissue' not in index:
                index.append('tissue')

            df_no_cat = self._junction[self._junction['delta_psi_cat'].isna()]
            self._junction_cat_concat = get_abs_max_rows(
                self._junction, index, 'delta_psi_cat') \
                .reset_index('event_type')
            self._junction_cat_concat = pd.concat(
                [self._junction_cat_concat, df_no_cat])

        return self._junction_cat_concat

    @property
    def splice_site_cat_concat(self):
        if self._splice_site_cat_concat is None:
            if self._splice_site is None:
                self.splice_site
            if 'tissue_cat' not in self._splice_site:
                raise IndexError(
                    '"tissue_cat not in columns. You have to run infer_cat first"')

            index = ['splice_site', 'event_type']
            if 'samples' in self.df:
                index.append('sample')
            if 'tissue' in self.df and 'tissue' not in index:
                index.append('tissue')

            df_no_cat = self._splice_site[self._splice_site['delta_psi_cat'].isna(
            )]
            self._splice_site_cat_concat = get_abs_max_rows(
                self._splice_site, index, 'delta_psi_cat') \
                .reset_index('event_type')
            self._splice_site_cat_concat = pd.concat(
                [self._splice_site_cat_concat, df_no_cat])

        return self._splice_site_cat_concat

    @property
    def gene_cat_concat(self):
        if self._gene_cat_concat is None:
            if self._gene is None:
                self.gene
            if 'tissue_cat' not in self._gene:
                raise IndexError(
                    '"tissue_cat not in columns. You have to run infer_cat first"')

            index = ['gene_name']
            if 'samples' in self.df:
                index.append('sample')
            if 'tissue' in self.df and 'tissue' not in index:
                index.append('tissue')

            df_no_cat = self._gene[self._gene['delta_psi_cat'].isna()]
            self._gene_cat_concat = get_abs_max_rows(
                self._gene[~self._gene['delta_psi_cat'].isna()], index, 'delta_psi_cat')
            self._gene_cat_concat = pd.concat(
                [self._gene_cat_concat, df_no_cat])

            if self.df_spliceAI is not None:
                self._gene_cat_concat = self._join_spliceAI(
                    self._gene_cat_concat)

        return self._gene_cat_concat

    @property
    def junction_cat_features(self):
        if self._junction_cat_features is None:
            if self._junction is None:
                self.junction
            if 'tissue_cat' not in self._junction:
                raise IndexError(
                    '"tissue_cat not in columns. You have to run infer_cat first"')

            cols_cat = [x for x in self._junction.columns if 'cat' in x]
            df_cat = self._junction[cols_cat]
            df_cat = df_cat.pivot(columns='tissue_cat')
            df_cat.columns = [f'{col[0]}'.replace('cat', f'{col[1]}') for col in df_cat.columns]

            df_target = self._junction[self._junction.columns.difference(
                cols_cat)]
            index = df_target.index.names
            df_target = df_target.reset_index().drop_duplicates(
            ).set_index(index)  # TODO: refactor this
            self._junction_cat_features = df_target.join(df_cat)

        return self._junction_cat_features

    @property
    def splice_site_cat_features(self):
        if self._splice_site_cat_features is None:
            if self._splice_site is None:
                self.splice_site
            if 'tissue_cat' not in self._splice_site:
                raise IndexError(
                    '"tissue_cat not in columns. You have to run infer_cat first"')

            cols_cat = [x for x in self._splice_site.columns if 'cat' in x]
            df_cat = self._splice_site[cols_cat]
            df_cat = df_cat.pivot(columns='tissue_cat')
            df_cat.columns = [f'{col[0]}'.replace('cat', f'{col[1]}') for col in df_cat.columns]

            df_target = self._splice_site[self._splice_site.columns.difference(
                cols_cat)]
            index = df_target.index.names
            df_target = df_target.reset_index().drop_duplicates(
            ).set_index(index)  # TODO: refactor this
            self._splice_site_cat_features = df_target.join(df_cat)

        return self._splice_site_cat_features

    @property
    def gene_cat_features(self):
        if self._gene_cat_features is None:
            if self._gene is None:
                self.gene
            if 'tissue_cat' not in self._gene:
                raise IndexError(
                    '"tissue_cat not in columns. You have to run infer_cat first"')

            cols_cat = [x for x in self._gene.columns if 'cat' in x]
            df_cat = self._gene[cols_cat]
            df_cat = df_cat.pivot(columns='tissue_cat')
            df_cat.columns = [f'{col[0]}'.replace('cat', f'{col[1]}') for col in df_cat.columns]

            df_target = self._gene[self._gene.columns.difference(cols_cat)]
            index = df_target.index.names
            df_target = df_target.reset_index().drop_duplicates(
            ).set_index(index)  # TODO: refactor this
            self._gene_cat_features = df_target.join(df_cat)

            if self.df_spliceAI is not None:
                self._gene_cat_features = self._join_spliceAI(
                    self._gene_cat_features)

        return self._gene_cat_features

    @staticmethod
    def _add_maf(df, population, default=-1):
        df['maf'] = df['variant'].map(
            lambda x: population.get(x, default))
        return df

    def add_maf(self, population, default=-1):
        self.df = self._add_maf(self.df, population, default)

    @staticmethod
    def _filter_private(df, max_num_sample=2):
        return df[df['samples'].str.split(';').map(set, na_action='ignore').map(list, na_action='ignore')
                  .map(len) <= max_num_sample]

    def _add_filter_maf(self, df, population=None,
                        maf_cutoff=0.001, default=-1):
        df = self._add_maf(df, population, default)
        return df[df['maf'] <= maf_cutoff]

    def filter_maf(self, max_num_sample=2, population=None,
                   maf_cutoff=0.001, default=-1):
        df = self.df
        df_spliceAI = self.df_spliceAI

        if max_num_sample:
            df = self._filter_private(df, max_num_sample)
            if df_spliceAI is not None:
                df_spliceAI = self._filter_private(df_spliceAI, max_num_sample)

        if population:
            df = self._add_filter_maf(df, population, maf_cutoff, default)
            if df_spliceAI is not None:
                df_spliceAI = self._add_filter_maf(
                    df_spliceAI, population, maf_cutoff, default)

        return SplicingOutlierResult(df, df_spliceAI)

    def predict_ensemble(self, pickle_file, df, features):
        self._ensemble = df.copy()
        model = pickle.load(open(pickle_file, 'rb'))
        X_test = self._ensemble[features].fillna(0)
        self._ensemble['ensemble_pred'] = model.predict_proba(X_test)[:, 1]
        return self._ensemble
