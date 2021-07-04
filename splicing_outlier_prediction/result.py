from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
from splicing_outlier_prediction.utils import get_abs_max_rows
from splicing_outlier_prediction.cat_dataloader import CatInference


class SplicingOutlierResult:

    def __init__(self, df, df_spliceAI=None):
        self.df = df
        self.df_spliceAI = df_spliceAI
        self._gene_spliceAI = None
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

    def add_samples(self, var_samples_df):
        '''
        var_samples_df: dataframe with variant and sample columns (e.g. output of 'to_sample_csv' from kipoiseq.extractors.vcf_query)
        '''
        var_samples_df['samples'] = var_samples_df \
            .groupby('variant')['sample'] \
            .transform(lambda x: ';'.join(x))
        var_samples_df = var_samples_df[['variant', 'samples']] \
            .drop_duplicates()
        self.df = self.df.set_index('variant') \
            .join(var_samples_df.set_index('variant')) \
            .reset_index()

    def filter_samples_with_RNA_seq(self, samples_for_tissue):
        """
        samples_for_tissue: Dict, keys: tissues, values: samples with RNA-seq for respective tissue

        filters out rows in self.df and self.df_spliceAI which do not have RNA-seq information
        """
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
        """
        max aggregate 'delta_psi' scores (MMSplice with SpliceMap) on junction level (if provided, for each sample and tissue)
        stores scores in self._junction
        """
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
        """
        max aggregate 'delta_psi' scores (MMSplice with SpliceMap) on splice site level (if provided, for each sample and tissue)
        stores scores in self._splice_site
        """
        if self._splice_site is None:
            index = ['splice_site', 'event_type']
            if 'samples' in self.df:
                index.append('sample')
            if 'tissue' in self.df:
                index.append('tissue')
            self._splice_site = get_abs_max_rows(
                self.junction, index, 'delta_psi') \
                .reset_index('event_type')
        return self._splice_site

    @property
    def gene(self):
        """
        max aggregate 'delta_psi' scores (MMSplice with SpliceMap) on gene level (if provided, for each sample and tissue)
        stores scores in self._gene
        """
        if self._gene is None:
            index = ['gene_name']
            if 'samples' in self.df:
                index.append('sample')
            if 'tissue' in self.df and 'tissue' not in index:
                index.append('tissue')
            self._gene = get_abs_max_rows(
                self.junction, index, 'delta_psi')

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
        """
        checks if self.df_spliceAI has 'tissue' column.
        If self.df has 'tissue' column and self.df_spliceAI does not have 'tissue' column,
        tissue independent spliceAI predictions are copied for each tissue in self.df
        """
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

    @property
    def gene_spliceAI(self):
        """
        max aggregate spliceAI scores on gene level (if provided, for each sample and tissue)
        stores scores in self._gene_spliceAI
        """
        # TODO: sample handling of SpliceAI is not good.
        if self._gene_spliceAI is None:
            index_spliceAI = ['gene_name']
            if 'tissue' in self.df.columns and 'tissue' not in self.df_spliceAI.columns:
                df_spliceAI, index_spliceAI = self._add_tissue_info_to_spliceAI()
            else:
                df_spliceAI = self.df_spliceAI
                if 'tissue' in self.df.columns:
                    index_spliceAI.append('tissue')
            if 'samples' in self.df_spliceAI:
                index_spliceAI.append('sample')
                df_spliceAI = self._explode_samples(df_spliceAI)
            index_spliceAI = sorted(index_spliceAI)
            self._gene_spliceAI = get_abs_max_rows(
                df_spliceAI, index_spliceAI, 'delta_score')
        return self._gene_spliceAI

    def _join_spliceAI_gene(self, df):
        self._gene_spliceAI = self.gene_spliceAI
        return df.join(self._gene_spliceAI, how='outer', rsuffix='_spliceAI')

    def infer_cat(self, cat_inference, progress=False):
        """
        cat_inference: List[CatInference] or CatInference

        infers delta_score_cat for each cat tissue in each target tissue, 
        based on ref_psi_target and measured delta_logit_psi in cat
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
                if cat.contains(junction, sample, tissue, row['event_type']):
                    infer_rows.append(
                        cat.infer(junction, sample, tissue, row['event_type']))
        df = pd.DataFrame(infer_rows)
        df = df.drop_duplicates().set_index(['junction', 'sample', 'tissue'])
        # self._junction can contain multiple cats (junction, sample, tissue) is not unique index
        self._junction = self.junction.join(df)
        self._splice_site = None
        self._gene = None

    def _cat_max_agg_level(self, index):
        """
        index: groupby for get_abs_max_rows (e.g. index=['gene_name', 'tissue'])

        function performs max aggregation independently for DNA based predictions ('delta_psi')
        and inferred scores from CAT ('delta_psi_cat'). Max aggregated scores are then joined.
        (Independent max aggregation is performed because the junction with maximum 'delta_psi' 
        could have missing value for 'delta_psi_cat' even though there is another junction 
        (on the same level, e.g. same gene) with large 'delta_psi_cat', but slightly smaller 'delta_psi'))
        """
        if self._junction is None or 'tissue_cat' not in self._junction:
            raise IndexError('"You have to run infer_cat first"')
        if 'samples' in self.df:
            index.append('sample')
        if 'tissue' in self.df:
            index.append('tissue')
        df_dna = get_abs_max_rows(
            self._junction, index, 'delta_psi')
        df_cat = get_abs_max_rows(
            self._junction[~self._junction['delta_psi_cat'].isna()],
            index, 'delta_psi_cat')
        dna_cols = [col for col in df_dna.columns if 'cat' not in col]
        cat_cols = [col for col in df_dna.columns if 'cat' in col]
        return df_dna[dna_cols].join(df_cat[cat_cols])

    def _cat_features_pivot_level(self, index):
        """
        index: groupby for get_abs_max_rows (e.g. index=['gene_name', 'tissue'])

        function runs 'self._cat_max_agg_level' for each 'tissue_cat'. 
        Then pivot cat information; allows to use different 'tissue_cat' as independent features.
        """
        if 'tissue_cat' not in index:
            raise IndexError(
                '"tissue_cat has to be in index to pivot dataframe"')
        self._junction['tissue_cat'] = self._junction['tissue_cat'].astype(str)
        df = self._cat_max_agg_level(index)
        self._junction['tissue_cat'] = self._junction['tissue_cat'].replace(
            'nan', np.nan)
        df = df.reset_index('tissue_cat')

        cols_cat = [x for x in df.columns if 'cat' in x]
        df_cat = df[cols_cat]
        df_cat = df_cat.pivot(columns='tissue_cat')
        df_cat.columns = [f'{col[0]}'.replace('cat', f'{col[1]}') for col in df_cat.columns]

        df_target = df[df.columns.difference(cols_cat)]
        index.remove('tissue_cat')
        df_target = df_target.reset_index().drop_duplicates(
        ).set_index(index)  # TODO: refactor this
        return df_target.join(df_cat)

    @property
    def junction_cat_concat(self):
        """
        run 'self._cat_max_agg_level' on junction level
        """
        if self._junction_cat_concat is None:
            index = ['junction', 'event_type']
            self._junction_cat_concat = self._cat_max_agg_level(
                index).reset_index('event_type')
        return self._junction_cat_concat

    @property
    def splice_site_cat_concat(self):
        """
        run 'self._cat_max_agg_level' on splice site level
        """
        if self._splice_site_cat_concat is None:
            index = ['splice_site', 'event_type']
            self._splice_site_cat_concat = self._cat_max_agg_level(
                index).reset_index('event_type')
        return self._splice_site_cat_concat

    @property
    def gene_cat_concat(self):
        """
        run 'self._cat_max_agg_level' on gene level
        Additionally join with max aggregated spliceAI scores (if provided)
        """
        if self._gene_cat_concat is None:
            index = ['gene_name']
            self._gene_cat_concat = self._cat_max_agg_level(index)

            if self.df_spliceAI is not None:
                self._gene_cat_concat = self._join_spliceAI_gene(
                    self._gene_cat_concat)

        return self._gene_cat_concat

    @property
    def junction_cat_features(self):
        """
        run 'self._cat_features_pivot_level' on junction level
        """
        if self._junction_cat_features is None:
            index = ['junction', 'event_type', 'tissue_cat']
            self._junction_cat_features = self._cat_features_pivot_level(
                index).reset_index('event_type')
        return self._junction_cat_features

    @property
    def splice_site_cat_features(self):
        """
        run 'self._cat_features_pivot_level' on splice site level
        """
        if self._splice_site_cat_features is None:
            index = ['splice_site', 'event_type', 'tissue_cat']
            self._splice_site_cat_features = self._cat_features_pivot_level(
                index).reset_index('event_type')
        return self._splice_site_cat_features

    @property
    def gene_cat_features(self):
        """
        run 'self._cat_features_pivot_level' on gene level
        Additionally join with max aggregated spliceAI scores (if provided)
        """
        if self._gene_cat_features is None:
            index = ['gene_name', 'tissue_cat']
            self._gene_cat_features = self._cat_features_pivot_level(index)

            if self.df_spliceAI is not None:
                self._gene_cat_features = self._join_spliceAI_gene(
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
        """
        pickle_file: path to pickle file with pre trained model
        df: input dataframe (contains features to be used for running model)
        features: features for running model (has to be in same order as in pre trained model)

        runs pre trained model (e.g. purely DNA based, or DNA+CAT) on features of SpliceOutlierResult
        """
        self._ensemble = df.copy()
        model = pickle.load(open(pickle_file, 'rb'))
        X_test = self._ensemble[features].fillna(0)
        self._ensemble['ensemble_pred'] = model.predict_proba(X_test)[:, 1]
        return self._ensemble
