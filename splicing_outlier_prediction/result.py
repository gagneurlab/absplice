from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
from pathlib import Path
import pathlib
from splicing_outlier_prediction.utils import get_abs_max_rows, normalize_gene_annotation, read_csv
from splicing_outlier_prediction.cat_dataloader import CatInference

class SplicingOutlierResult:

    def __init__(self, df_mmsplice=None, df_spliceai=None, df_mmsplice_cat=None, gene_map=None, gene_tpm=None):
        self.df_mmsplice = self.validate_df_mmsplice(df_mmsplice)
        self.df_mmsplice_cat = self.validate_df_mmsplice_cat(df_mmsplice_cat)
        self.gene_map = self.validate_df_gene_map(gene_map)
        self.gene_tpm = self.validate_df_gene_tpm(gene_tpm)
        self.df_spliceai = self.validate_df_spliceai(df_spliceai)
        self._df_spliceai_tissue = None
        self.contains_chr = self._contains_chr()
        self._junction = None
        self._splice_site = None
        self._gene_mmsplice = None
        self._gene_mmsplice_cat = None
        self._gene_spliceai = None
        self._absplice_dna_input = None
        self._absplice_rna_input = None
        self._absplice_dna = None
        self._absplice_rna = None
        self._gene_absplice_dna = None
        self._gene_absplice_rna = None
    
    def _validate_df(self, df, columns):
        if not isinstance(df, pd.DataFrame):
            df = read_csv(df)
        assert pd.Series(columns).isin(df.columns).all()
        return df
    
    def validate_df_mmsplice(self, df_mmsplice):
        if df_mmsplice is not None:
            df_mmsplice = self._validate_df(
                df_mmsplice, 
                columns=['variant', 'gene_id', 'tissue', 
                         'delta_psi', 'ref_psi', 'median_n', 'gene_tpm'])
        return df_mmsplice
    
    def validate_df_mmsplice_cat(self, df_mmsplice_cat):
        if df_mmsplice_cat is not None:
            df_mmsplice_cat = self._validate_df(
                df_mmsplice_cat, 
                columns=['variant', 'gene_id', 'tissue', 
                         'delta_psi', 'ref_psi', 'median_n', 'gene_tpm', 
                         'tissue_cat', 'delta_psi_cat'])
        return df_mmsplice_cat
    
    def validate_df_spliceai(self, df_spliceai):
        if df_spliceai is not None:
            df_spliceai = self._validate_df(
                df_spliceai, 
                columns=['variant', 'gene_name', 'delta_score'])
            if self.gene_map is not None:
                df_spliceai = normalize_gene_annotation(df_spliceai, self.gene_map, key='gene_name', value='gene_id') 
        return df_spliceai
    
    def validate_df_gene_tpm(self, gene_tpm):
        if gene_tpm is not None:
            gene_tpm = self._validate_df(
                gene_tpm, 
                columns=['gene_id', 'tissue', 'gene_tpm'])
        if gene_tpm is not None and self.df_mmsplice is not None:
            missing_tissues = set(self.df_mmsplice['tissue']).difference(set(gene_tpm['tissue']))
            if len(missing_tissues) > 0:
                raise KeyError(" %s are missing in gene_tpm" % missing_tissues)
        return gene_tpm
    
    def validate_df_gene_map(self, gene_map):
        if gene_map is not None:
            gene_map = self._validate_df(
                gene_map, 
                columns=['gene_id', 'gene_name'])
        return gene_map
        
    def _contains_chr(self):
        if self.df_mmsplice is not None:
            return 'chr' in self.df_mmsplice.junction[0]
        else:
            return None
        
    def add_spliceai(self, df, 
                     gene_mapping=True, key='gene_name', value='gene_id'):
        """
        Includes spliceai predictions into results.

        Args:
          df: path to csv or dataframe of spliceai predictions
        """
        if type(df) == str:
            df = pd.read_csv(df)
        self.df_spliceai = df
        if gene_mapping:
            self.df_spliceai = normalize_gene_annotation(self.df_spliceai, self.gene_map, key, value)
    
    def _add_tissue_info_to_spliceai(self):
        """
        checks if self.df_spliceai has 'tissue' column.
        If self.df_mmsplice has 'tissue' column and self.df_spliceai does not have 'tissue' column,
        tissue independent spliceai predictions are copied for each tissue in self.df_mmsplice
        """
        df_spliceai = self.df_spliceai
        if self.gene_tpm is not None: #only duplicate predictions for expressed genes
            df_gene_tpm =  self.gene_tpm[
                self.gene_tpm['tissue'].isin(self.df_mmsplice['tissue'].unique())]
            self._df_spliceai_tissue = df_spliceai.set_index('gene_id')\
                .join(df_gene_tpm.set_index('gene_id'))\
                .reset_index()
        else: #if gene tpm values not provided for tissues, duplicate all predictions
            l = list()
            for tissue in self.df_mmsplice['tissue'].unique():
                _df = df_spliceai.copy()
                _df['tissue'] = tissue
                l.append(_df)
            self._df_spliceai_tissue = pd.concat(l)
        return self._df_spliceai_tissue

    def _add_samples(self, df, var_samples_df):
        var_samples_df = var_samples_df[['variant', 'sample']] \
            .drop_duplicates()
        df = df.set_index('variant') \
            .join(var_samples_df.set_index('variant'),
                  how='inner') \
            .reset_index()
        return df

    def add_samples(self, var_samples_df):
        '''
        var_samples_df: dataframe with variant and sample columns 
        (e.g. output of 'to_sample_csv' from kipoiseq.extractors.vcf_query)
        '''
        if self.df_mmsplice is not None:
            if 'sample' not in self.df_mmsplice.columns:
                self.df_mmsplice = self._add_samples(self.df_mmsplice, var_samples_df)
        if self.df_spliceai is not None:
            if 'sample' not in self.df_spliceai.columns:
                self.df_spliceai = self._add_samples(self.df_spliceai, var_samples_df)
            
    @staticmethod
    def _explode(df, col='samples', new_name=None):
        if new_name == None:
            new_name = col
        df = df.copy()
        df[col] = df[col].str.split(';').map(
            set, na_action='ignore').map(list, na_action='ignore')
        return df.rename(columns={col: new_name}).explode(new_name)
    
    # def _filter_samples_with_RNA_seq(self, df, samples_for_tissue):
    #     l = list()
    #     for tissue, samples in samples_for_tissue.items():
    #         df_with_RNA = df[(df['tissue'] == tissue) & (df['sample'].isin(samples))]
    #         l.append(df_with_RNA)
    #     return pd.concat(l)

    # def filter_samples_with_RNA_seq(self, samples_for_tissue):
    #     """
    #     samples_for_tissue: Dict, keys: tissue, values: samples with RNA-seq for respective tissue

    #     filters out rows in self.df_mmsplice and self.df_spliceai which do not have RNA-seq information
    #     """
    #     if self.df_mmsplice is not None:
    #         if 'tissue' not in self.df_mmsplice or 'sample' not in self.df_mmsplice:
    #             raise KeyError('tissue or sample annotation missing')
    #         self.df_mmsplice = self._filter_samples_with_RNA_seq(self.df_mmsplice, samples_for_tissue)

    #     if self.df_spliceai is not None:
    #         if 'tissue' not in self.df_spliceai:
    #             self.df_spliceai = self._add_tissue_info_to_spliceai()
    #         self.df_spliceai = self._filter_samples_with_RNA_seq(self.df_spliceai, samples_for_tissue)

    def infer_cat(self, cat_inference, progress=False):
        """
        cat_inference: List[CatInference] or CatInference

        infers delta_score_cat for each cat tissue in each target tissue, 
        based on ref_psi_target and measured delta_logit_psi in cat
        """
        if 'sample' not in self.df_mmsplice.columns:
            raise ValueError('"sample" column is missing. Call add.samples() first')

        if type(cat_inference) == CatInference:
            cat_inference = [cat_inference]

        infer_rows = list()
        for cat in cat_inference:
            assert self.contains_chr == cat.contains_chr
            rows = self.junction.iterrows()
            if progress:
                rows = tqdm(rows, total=self.junction.shape[0])
            for (junction, tissue, sample), row in rows:
                if cat.contains(junction, tissue, sample, row['event_type']):
                    infer_rows.append(
                        cat.infer(junction, tissue, sample, row['event_type']))
        df = pd.DataFrame(infer_rows)
        assert df.shape[0] > 0
        df = df.drop_duplicates().set_index(['junction', 'tissue', 'sample'])
        # self._junction can contain multiple cats (junction, tissue, sample) is not unique index
        self._junction = self.junction.join(df)
        self.df_mmsplice_cat = self._junction[~self._junction['tissue_cat'].isna()]
        # self._splice_site = None
        # self._gene = None
        
    def _get_maximum_effect(self, df, groupby, score, dropna=True):
        if not isinstance(df.index, pd.RangeIndex):
            df = df.reset_index()
        assert len(set(groupby).difference(df.columns)) == 0
        return get_abs_max_rows(df.set_index(groupby), groupby, score, dropna)

    @property
    def psi5(self):
        return SplicingOutlierResult(self.df_mmsplice[self.df_mmsplice['event_type'] == 'psi5'])

    @property
    def psi3(self):
        return SplicingOutlierResult(self.df_mmsplice[self.df_mmsplice['event_type'] == 'psi3'])
    
    @property
    def junction(self): #NOTE: max aggregate over all variants
        groupby=['junction', 'event_type', 'tissue']
        if 'sample' in self.df_mmsplice:
            groupby.append('sample')
        if self._junction is None:
            self._junction = self._get_maximum_effect(self.df_mmsplice, groupby, score='delta_psi')
            self._junction = self._junction.reset_index('event_type')
        return self._junction

    @property
    def splice_site(self): #NOTE: max aggregate over all variants
        groupby=['splice_site', 'event_type', 'tissue']
        if 'sample' in self.df_mmsplice:
            groupby.append('sample')
        if self._splice_site is None:
            self._splice_site = self._get_maximum_effect(self.df_mmsplice, groupby, score='delta_psi')
            self._splice_site = self._splice_site.reset_index('event_type')
        return self._splice_site

    @property
    def gene_mmsplice(self): #NOTE: max aggregate over all variants
        groupby=['gene_id', 'tissue']
        if 'sample' in self.df_mmsplice:
            groupby.append('sample')
        if self._gene_mmsplice is None:
            # TODO: if cat infer was called, add it to gene level aggregation
            self._gene_mmsplice = self._get_maximum_effect(self.df_mmsplice, groupby, score='delta_psi')
        return self._gene_mmsplice
    
    @property
    def gene_mmsplice_cat(self): #NOTE: max aggregate over all variants
        groupby=['gene_id', 'tissue']
        if 'sample' in self.df_mmsplice_cat:
            groupby.append('sample')
        if self._gene_mmsplice_cat is None:
            self._gene_mmsplice_cat = self._get_maximum_effect(self.df_mmsplice_cat, groupby, score='delta_psi_cat')
        return self._gene_mmsplice_cat

    @property
    def gene_spliceai(self): #NOTE: max aggregate over all variants
        groupby=['gene_id']
        if 'sample' in self.df_spliceai:
            groupby.append('sample')
        if self._gene_spliceai is None:
            self._gene_spliceai = self._get_maximum_effect(self.df_spliceai, groupby, score='delta_score')
        return self._gene_spliceai

    @property
    def absplice_dna_input(self):
        groupby=['variant', 'gene_id', 'tissue']
        if 'sample' in self.df_mmsplice and 'sample' in self.df_spliceai:
            groupby.append('sample')
        if self._absplice_dna_input is None: 
            df_spliceai = self._add_tissue_info_to_spliceai() # add tissue info to spliceai
            df_mmsplice = self._get_maximum_effect(self.df_mmsplice, groupby, score='delta_psi')
            df_spliceai = self._get_maximum_effect(df_spliceai, groupby, score='delta_score', dropna=False) #dropna=False assures that also missing gene_id and genes that do not have tpm values in tissues are predicted
            cols_spliceai = ['delta_score', 'gene_name']
            cols_mmsplice = [
                'Chromosome', 'Start', 'End', 'Strand', 'junction', 'event_type', 'splice_site',
                'delta_logit_psi', 'delta_psi', 'ref_psi', 'k', 'n', 'median_n', 'gene_tpm', 'gene_name',
                'novel_junction', 'weak_site_donor', 'weak_site_acceptor']
            self._absplice_dna_input = df_mmsplice[cols_mmsplice].join(df_spliceai[cols_spliceai], how='outer', rsuffix='_spliceai')
        return self._absplice_dna_input
    
    @property
    def absplice_rna_input(self):
        groupby=['variant', 'gene_id', 'tissue']
        if 'sample' in self.df_mmsplice and 'sample' in self.df_spliceai:
            groupby.append('sample')
        if self._absplice_rna_input is None: 
            df_mmsplice_cat = self._get_maximum_effect(self.df_mmsplice_cat, groupby, score='delta_psi_cat') 
            cols_mmsplice_cat = [
                'junction', 'delta_psi', 'ref_psi', 'median_n', 
                *[col for col in df_mmsplice_cat.columns if 'cat' in col]]
            self._absplice_rna_input = self.absplice_dna_input.join(df_mmsplice_cat[cols_mmsplice_cat], how='outer', rsuffix='_from_cat_infer')
        return self._absplice_rna_input
    
    def _predict_absplice(self, df, absplice_score, pickle_file, features, abs_features, median_n_cutoff, tpm_cutoff):
        model = pickle.load(open(pickle_file, 'rb'))
        df['splice_site_is_expressed'] = df['median_n'] > median_n_cutoff
        df['gene_is_expressed'] = df['gene_tpm'] > tpm_cutoff
        df = df[features].fillna(0)
        if abs_features == True:
            df = np.abs(df)
        df[absplice_score] = model.predict_proba(df)[:, 1]
        return df
        
    def predict_absplice_dna(self, pickle_file=None, abs_features=True, median_n_cutoff=0, tpm_cutoff=1):
        features = ['delta_psi', 'delta_score',
                    'gene_is_expressed', 'splice_site_is_expressed',
                    'ref_psi', 'delta_logit_psi']
        if pickle_file == None:
            pickle_file = pickle_absplice_DNA
            
        self.absplice_dna = self._predict_absplice(
            df=self.absplice_dna_input,
            absplice_score='AbSplice_DNA', 
            pickle_file=pickle_file, 
            features=features, 
            abs_features=abs_features, 
            median_n_cutoff=median_n_cutoff, 
            tpm_cutoff=tpm_cutoff)
        return self.absplice_dna
    
    def predict_absplice_rna(self, pickle_file=None, abs_features=True, median_n_cutoff=0, tpm_cutoff=1):
        features = ['delta_psi', 'delta_score', 
                    'gene_is_expressed', 'splice_site_is_expressed', 
                    'ref_psi', 'delta_logit_psi',
                    'delta_psi_cat', 'count_cat',
                    'psi_cat', 'ref_psi_cat']
        if pickle_file == None:
            pickle_file = pickle_absplice_RNA
            
        self.absplice_rna = self._predict_absplice(
            df=self.absplice_rna_input,
            absplice_score='AbSplice_RNA', 
            pickle_file=pickle_file, 
            features=features, 
            abs_features=abs_features, 
            median_n_cutoff=median_n_cutoff, 
            tpm_cutoff=tpm_cutoff)
        return self.absplice_rna
    
    @property
    def gene_absplice_dna(self): #NOTE: max aggregate over all variants
        groupby=['gene_id', 'tissue']
        if 'sample' in self.df_mmsplice and 'sample' in self.df_spliceai:
            groupby.append('sample')
        if self._gene_absplice_dna is None:
            self._gene_absplice_dna = self._get_maximum_effect(self.absplice_dna, groupby, score='AbSplice_DNA')
        return self._gene_absplice_dna
    
    @property
    def gene_absplice_rna(self): #NOTE: max aggregate over all variants
        groupby=['gene_id', 'tissue']
        if 'sample' in self.df_mmsplice and 'sample' in self.df_spliceai:
            groupby.append('sample')
        if self._gene_absplice_rna is None:
            self._gene_absplice_rna = self._get_maximum_effect(self.absplice_rna, groupby, score='AbSplice_RNA')
        return self._gene_absplice_rna
    
    @staticmethod
    def _add_maf(df, population, default=-1):
        df['maf'] = df['variant'].map(
            lambda x: population.get(x, default))
        return df

    @staticmethod
    def _filter_private(df, max_num_sample=2):
        df['samples'] = df.groupby('variant')['sample'].apply(lambda x: ';'.join(x))
        return df[df['samples'].str.split(';').map(set, na_action='ignore').map(list, na_action='ignore')
                  .map(len) <= max_num_sample]

    def _add_filter_maf(self, df, population=None,
                        maf_cutoff=0.001, default=-1):
        df = self._add_maf(df, population, default)
        return df[df['maf'] <= maf_cutoff]

    def filter_maf(self, max_num_sample=2, population=None,
                   maf_cutoff=0.001, default=-1):
        df_mmsplice = self.df_mmsplice
        df_spliceai = self.df_spliceai

        if max_num_sample:
            df_mmsplice = self._filter_private(df_mmsplice, max_num_sample)
            if df_spliceai is not None:
                df_spliceai = self._filter_private(df_spliceai, max_num_sample)

        if population:
            df_mmsplice = self._add_filter_maf(df_mmsplice, population, maf_cutoff, default)
            if df_spliceai is not None:
                df_spliceai = self._add_filter_maf(
                    df_spliceai, population, maf_cutoff, default)

        return SplicingOutlierResult(df_mmsplice, df_spliceai)

    

###### FROM BEFORE

# @property
# def gene_rna(self):
#     groupby=['variant', 'gene_id', 'tissue']
#     if 'sample' in self.df_mmsplice and 'sample' in self.df_spliceai:
#         groupby.append('sample')
#     if self._gene_rna is None:
#         # get maximum variant effect on gene
#         if self.df_mmsplice_cat:
#             df_mmsplice_cat = self._get_maximum_effect(self.df_mmsplice_cat, groupby, score='delta_psi_cat') # maybe there has to be tissue_cat in index as well (otherwise max agg over all tissue_cats)
#         else:
#             raise AttributeError("You have to call infer_cat first")
#         # dna_cols = [col for col in self._absplice_input.columns if 'cat' not in col]
#         cat_cols = [col for col in df_mmsplice_cat.columns if 'cat' in col]
#         self._gene_rna = self._absplice_input.join(df_mmsplice_cat[cat_cols], how='outer', rsuffix='_from_cat_infer')
#     return self._gene_rna

# def _join_spliceai_gene(self, df):
#     self._gene_spliceai = self.gene_spliceai
#     return df.join(self._gene_spliceai, how='outer', rsuffix='_spliceai')
    
    
# def _cat_max_agg_level(self, index):
#     """
#     index: groupby for get_abs_max_rows (e.g. index=['gene_name', 'tissue'])

#     function performs max aggregation independently for DNA based predictions ('delta_psi')
#     and inferred scores from CAT ('delta_psi_cat'). Max aggregated scores are then joined.
#     (Independent max aggregation is performed because the junction with maximum 'delta_psi' 
#     could have missing value for 'delta_psi_cat' even though there is another junction 
#     (on the same level, e.g. same gene) with large 'delta_psi_cat', but slightly smaller 'delta_psi'))
#     """
#     if self._junction is None or 'tissue_cat' not in self._junction:
#         raise IndexError('"You have to run infer_cat first"')
#     if 'sample' in self.df_mmsplice:
#         index.append('sample')
#     if 'tissue' in self.df_mmsplice:
#         index.append('tissue')
#     df_dna = get_abs_max_rows(
#         self._junction, index, 'delta_psi')
#     df_cat = get_abs_max_rows(
#         self._junction[~self._junction['delta_psi_cat'].isna()],
#         index, 'delta_psi_cat')
#     dna_cols = [col for col in df_dna.columns if 'cat' not in col]
#     cat_cols = [col for col in df_dna.columns if 'cat' in col]
#     return df_dna[dna_cols].join(df_cat[cat_cols])

# def _cat_features_pivot_level(self, index):
#     """
#     index: groupby for get_abs_max_rows (e.g. index=['gene_name', 'tissue'])

#     function runs 'self._cat_max_agg_level' for each 'tissue_cat'. 
#     Then pivot cat information; allows to use different 'tissue_cat' as independent features.
#     """
#     if 'tissue_cat' not in index:
#         raise IndexError(
#             '"tissue_cat has to be in index to pivot dataframe"')
#     self._junction['tissue_cat'] = self._junction['tissue_cat'].astype(str)
#     df = self._cat_max_agg_level(index)
#     self._junction['tissue_cat'] = self._junction['tissue_cat'].replace(
#         'nan', np.nan)
#     df = df.reset_index('tissue_cat')

#     cols_cat = [x for x in df.columns if 'cat' in x]
#     df_cat = df[cols_cat]
#     df_cat = df_cat.pivot(columns='tissue_cat')
#     df_cat.columns = [f'{col[0]}'.replace('cat', f'{col[1]}') for col in df_cat.columns]

#     df_target = df[df.columns.difference(cols_cat)]
#     index.remove('tissue_cat')
#     df_target = df_target.reset_index().drop_duplicates(
#     ).set_index(index)  # TODO: refactor this
#     return df_target.join(df_cat)

# @property
# def junction_cat_concat(self):
#     """
#     run 'self._cat_max_agg_level' on junction level
#     """
#     if self._junction_cat_concat is None:
#         index = ['junction', 'event_type']
#         self._junction_cat_concat = self._cat_max_agg_level(
#             index).reset_index('event_type')
#     return self._junction_cat_concat

# @property
# def splice_site_cat_concat(self):
#     """
#     run 'self._cat_max_agg_level' on splice site level
#     """
#     if self._splice_site_cat_concat is None:
#         index = ['splice_site', 'event_type']
#         self._splice_site_cat_concat = self._cat_max_agg_level(
#             index).reset_index('event_type')
#     return self._splice_site_cat_concat

# @property
# def gene_cat_concat(self):
#     """
#     run 'self._cat_max_agg_level' on gene level
#     Additionally join with max aggregated spliceai scores (if provided)
#     """
#     if self._gene_cat_concat is None:
#         index = ['gene_name']
#         self._gene_cat_concat = self._cat_max_agg_level(index)

#         if self.df_spliceai is not None:
#             self._gene_cat_concat = self._join_spliceai_gene(
#                 self._gene_cat_concat)

#     return self._gene_cat_concat

# @property
# def junction_cat_features(self):
#     """
#     run 'self._cat_features_pivot_level' on junction level
#     """
#     if self._junction_cat_features is None:
#         index = ['junction', 'event_type', 'tissue_cat']
#         self._junction_cat_features = self._cat_features_pivot_level(
#             index).reset_index('event_type')
#     return self._junction_cat_features

# @property
# def splice_site_cat_features(self):
#     """
#     run 'self._cat_features_pivot_level' on splice site level
#     """
#     if self._splice_site_cat_features is None:
#         index = ['splice_site', 'event_type', 'tissue_cat']
#         self._splice_site_cat_features = self._cat_features_pivot_level(
#             index).reset_index('event_type')
#     return self._splice_site_cat_features

# @property
# def gene_cat_features(self):
#     """
#     run 'self._cat_features_pivot_level' on gene level
#     Additionally join with max aggregated spliceai scores (if provided)
#     """
#     if self._gene_cat_features is None:
#         index = ['gene_name', 'tissue_cat']
#         self._gene_cat_features = self._cat_features_pivot_level(index)

#         if self.df_spliceai is not None:
#             self._gene_cat_features = self._join_spliceai_gene(
#                 self._gene_cat_features)

#     return self._gene_cat_features