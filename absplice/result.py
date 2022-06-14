import resource
from pkg_resources import resource_filename
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
from pathlib import Path
import pathlib
from absplice.utils import get_abs_max_rows, normalize_gene_annotation, \
    read_csv, read_spliceai
from absplice.cat_dataloader import CatInference

GENE_MAP = resource_filename(
    'absplice', 'precomputed/GENE_MAP.tsv.gz')
GENE_TPM = resource_filename(
    'absplice', 'precomputed/GENE_TPM.csv.gz')
ABSPLICE_DNA = resource_filename(
    'absplice', 'precomputed/AbSplice_DNA.pkl')
ABSPLICE_RNA = resource_filename(
    'absplice', 'precomputed/AbSplice_RNA.pkl')

dtype_columns = {
    'variant': pd.StringDtype(),
    'gene_id': pd.StringDtype(),
    'tissue': pd.StringDtype(),
    'sample': pd.StringDtype(),
    'Chromosome': pd.StringDtype(),
    'Start': 'Int64',
    'End': 'Int64',
    'Strand': pd.StringDtype(),
    'junction': pd.StringDtype(),
    'event_type': pd.StringDtype(),
    'splice_site': pd.StringDtype(),
    'gene_name': pd.StringDtype(),
    'delta_logit_psi': 'float64',
    'delta_psi': 'float64',
    'ref_psi': 'float64',
    'k': 'Int64',
    'n': 'Int64',
    'median_n': 'float64',
    'novel_junction': pd.BooleanDtype(),
    'weak_site_donor': pd.BooleanDtype(),
    'weak_site_acceptor': pd.BooleanDtype(),
    'delta_score': 'float64',
    'gene_name': pd.StringDtype(),
    'gene_name_spliceai': pd.StringDtype(),
    'gene_tpm': 'float64',
    'tissue_cat':  pd.StringDtype(),
    'k_cat': 'Int64',
    'n_cat': 'Int64',
    'median_n_cat': 'float64',
    'psi_cat': 'float64',
    'ref_psi_cat': 'float64',
    'delta_logit_psi_cat': 'float64',
    'delta_psi_cat': 'float64',
    'AbSplice_DNA': 'float64',
    'AbSplice_RNA': 'float64',
}


class SplicingOutlierResult:

    def __init__(self,
                 df_mmsplice=None,
                 df_spliceai=None,
                 df_mmsplice_cat=None,
                 gene_map=None,
                 gene_tpm=None,
                 df_var_samples=None,
                 df_absplice_dna_input=None,
                 df_absplice_rna_input=None,
                 df_absplice_dna=None,
                 df_absplice_rna=None,
                 ):
        self.df_var_samples = self.validate_df_var_samples(df_var_samples)
        self.df_mmsplice = self.validate_df_mmsplice(df_mmsplice)
        self.df_mmsplice_cat = self.validate_df_mmsplice_cat(df_mmsplice_cat)
        self.gene_map = self.validate_df_gene_map(gene_map)
        self.gene_tpm = self.validate_df_gene_tpm(gene_tpm)
        self.df_spliceai = self.validate_df_spliceai(df_spliceai)
        self._absplice_dna_input = self.validate_absplice_dna_input(
            df_absplice_dna_input)
        self._absplice_rna_input = self.validate_absplice_rna_input(
            df_absplice_rna_input)
        self._absplice_dna = self.validate_absplice_dna(df_absplice_dna)
        self._absplice_rna = self.validate_absplice_rna(df_absplice_rna)
        self.contains_chr = self._contains_chr()
        self.contains_samples = self._contains_samples()
        self._df_spliceai_tissue = None
        self._junction = None
        self._splice_site = None
        self._gene_mmsplice = None
        self._gene_mmsplice_cat = None
        self._gene_spliceai = None
        self._gene_absplice_dna = None
        self._gene_absplice_rna = None
        self._variant_mmsplice = None
        self._variant_mmsplice_cat = None
        self._variant_spliceai = None
        self._variant_absplice_dna = None
        self._variant_absplice_rna = None

    def _validate_df(self, df, columns):
        if not isinstance(df, pd.DataFrame):
            df = read_csv(df)
        df = df.reset_index()
        if 'index' in df.columns:
            df = df.drop(columns='index')
        assert pd.Series(columns).isin(df.columns).all()
        return df

    def _validate_dtype(self, df):
        for col in df.columns:
            if col in dtype_columns.keys():
                df = df.astype({col: dtype_columns[col]})
        return df

    def validate_df_mmsplice(self, df_mmsplice):
        if df_mmsplice is not None:
            df_mmsplice = self._validate_df(
                df_mmsplice,
                columns=['variant', 'gene_id', 'tissue',
                         'delta_psi', 'ref_psi', 'median_n', 'gene_tpm'])
            df_mmsplice = self._validate_dtype(df_mmsplice)
            if self.df_var_samples is not None:
                df_mmsplice = self._add_samples(df_mmsplice)
        return df_mmsplice

    def validate_df_mmsplice_cat(self, df_mmsplice_cat):
        if df_mmsplice_cat is not None:
            df_mmsplice_cat = self._validate_df(
                df_mmsplice_cat,
                columns=['variant', 'gene_id', 'tissue',
                         'delta_psi', 'ref_psi', 'median_n', 'gene_tpm',
                         'tissue_cat', 'delta_psi_cat'])
            df_mmsplice_cat = self._validate_dtype(df_mmsplice_cat)
            df_mmsplice_cat = df_mmsplice_cat[
                ~df_mmsplice_cat['delta_psi_cat'].isna()
            ]
        return df_mmsplice_cat

    def validate_df_spliceai(self, df_spliceai):
        if df_spliceai is not None:
            df_spliceai = read_spliceai(df_spliceai)
            df_spliceai = self._validate_df(
                df_spliceai,
                columns=['variant', 'gene_name', 'delta_score'])
            if self.gene_map is not None:
                df_spliceai = normalize_gene_annotation(
                    df_spliceai, self.gene_map, key='gene_name', value='gene_id')
            df_spliceai = self._validate_dtype(df_spliceai)
            if self.df_var_samples is not None:
                df_spliceai = self._add_samples(df_spliceai)
        return df_spliceai

    def validate_df_var_samples(self, df_var_samples):
        if df_var_samples is not None:
            df_var_samples = self._validate_df(
                df_var_samples,
                columns=['variant', 'sample'])
            df_var_samples = self._validate_dtype(df_var_samples)
            df_var_samples = df_var_samples[[
                'variant', 'sample']].drop_duplicates()
        return df_var_samples

    def validate_df_gene_tpm(self, gene_tpm):
        if gene_tpm is not None:
            gene_tpm = self._validate_df(
                gene_tpm,
                columns=['gene_id', 'tissue', 'gene_tpm'])
            gene_tpm = gene_tpm[['gene_id', 'tissue', 'gene_tpm']]
            gene_tpm = self._validate_dtype(gene_tpm)
        else:
            gene_tpm = self._validate_df(
                GENE_TPM,
                columns=['gene_id', 'tissue', 'gene_tpm'])
        if gene_tpm is not None and self.df_mmsplice is not None:
            missing_tissues = set(self.df_mmsplice['tissue']).difference(
                set(gene_tpm['tissue']))
            if len(missing_tissues) > 0:
                raise KeyError(" %s are missing in gene_tpm" % missing_tissues)
        return gene_tpm

    def validate_df_gene_map(self, gene_map):
        if gene_map is not None:
            gene_map = self._validate_df(
                gene_map,
                columns=['gene_id', 'gene_name'])
            gene_map = self._validate_dtype(gene_map)
        else:
            gene_map = self._validate_df(
                GENE_MAP,
                columns=['gene_id', 'gene_name'])
        return gene_map

    def validate_absplice_dna_input(self, df_absplice_dna_input):
        if df_absplice_dna_input is not None:
            df_absplice_dna_input = self._validate_df(
                df_absplice_dna_input,
                columns=[
                    'variant', 'gene_id', 'tissue', 'gene_name', 'gene_name_spliceai', 'gene_tpm',
                    'Chromosome', 'Start', 'End', 'Strand', 'junction', 'event_type', 'splice_site',
                    'delta_score', 'delta_logit_psi', 'delta_psi', 'ref_psi', 'k', 'n', 'median_n',
                    'novel_junction', 'weak_site_donor', 'weak_site_acceptor'
                ])
            df_absplice_dna_input = self._validate_dtype(df_absplice_dna_input)
            groupby = ['variant', 'gene_id', 'tissue']
            if 'sample' in df_absplice_dna_input:
                groupby.append('sample')
            df_absplice_dna_input = df_absplice_dna_input.set_index(groupby)
        return df_absplice_dna_input

    def validate_absplice_rna_input(self, df_absplice_rna_input):
        if df_absplice_rna_input is not None:
            df_absplice_rna_input = self._validate_df(
                df_absplice_rna_input,
                columns=[
                    'variant', 'gene_id', 'tissue', 'sample', 'gene_name', 'gene_name_spliceai', 'gene_tpm',
                    'Chromosome', 'Start', 'End', 'Strand', 'junction', 'event_type', 'splice_site',
                    'delta_score', 'delta_logit_psi', 'delta_psi', 'ref_psi', 'k', 'n', 'median_n',
                    'novel_junction', 'weak_site_donor', 'weak_site_acceptor',
                    'tissue_cat', 'k_cat', 'n_cat', 'median_n_cat', 'psi_cat', 'ref_psi_cat',
                    'delta_logit_psi_cat', 'delta_psi_cat'
                ])
            df_absplice_rna_input = self._validate_dtype(df_absplice_rna_input)
            groupby = ['variant', 'gene_id', 'tissue', 'sample']
            df_absplice_rna_input = df_absplice_rna_input.set_index(groupby)
        return df_absplice_rna_input

    def validate_absplice_dna(self, df_absplice_dna):
        if df_absplice_dna is not None:
            df_absplice_dna = self._validate_df(
                df_absplice_dna,
                columns=[
                    'gene_id', 'tissue', 'AbSplice_DNA'
                ])
            df_absplice_dna = self._validate_dtype(df_absplice_dna)
        return df_absplice_dna

    def validate_absplice_rna(self, df_absplice_rna):
        if df_absplice_rna is not None:
            df_absplice_rna = self._validate_df(
                df_absplice_rna,
                columns=[
                    'gene_id', 'tissue', 'AbSplice_RNA'
                ])
            df_absplice_rna = self._validate_dtype(df_absplice_rna)
        return df_absplice_rna

    def _contains_chr(self):
        if self.df_mmsplice is not None:
            return 'chr' in self.df_mmsplice.junction[0]
        else:
            return None
        
    def _contains_samples(self):
        contains_samples = False
        if self.df_mmsplice is not None:
            if 'sample' in self.df_mmsplice.columns:
                contains_samples = True
        elif self.df_spliceai is not None:
            if 'sample' in self.df_spliceai.columns:
                contains_samples = True
        return contains_samples
        
    def _add_tissue_info_to_spliceai(self):
        """
        checks if self.df_spliceai has 'tissue' column.
        If self.df_mmsplice has 'tissue' column and self.df_spliceai does not have 'tissue' column,
        tissue independent spliceai predictions are copied for each tissue in self.df_mmsplice
        """
        df_spliceai = self.df_spliceai
        if self.df_mmsplice is not None:
            l = list()
            for tissue in self.df_mmsplice['tissue'].unique():
                _df = df_spliceai.copy()
                _df['tissue'] = tissue
                l.append(_df)
            self._df_spliceai_tissue = pd.concat(l)
        else:
            self._df_spliceai_tissue = df_spliceai.copy()
            self._df_spliceai_tissue['tissue'] = 'Not provided'
        return self._df_spliceai_tissue

    def _add_samples(self, df):
        df = df.set_index('variant') \
            .join(self.df_var_samples.set_index('variant'),
                  how='inner') \
            .reset_index()
        return df

    def add_samples(self, df_var_samples):
        '''
        df_var_samples: dataframe with variant and sample columns 
        (e.g. output of 'to_sample_csv' from kipoiseq.extractors.vcf_query)
        '''
        self.df_var_samples = self.validate_df_var_samples(df_var_samples)
        if self.df_mmsplice is not None:
            if 'sample' not in self.df_mmsplice.columns:
                self.df_mmsplice = self._add_samples(self.df_mmsplice)
        if self.df_spliceai is not None:
            if 'sample' not in self.df_spliceai.columns:
                self.df_spliceai = self._add_samples(self.df_spliceai)

    def infer_cat(self, cat_inference, progress=False):
        """
        cat_inference: List[CatInference] or CatInference

        infers delta_score_cat for each cat tissue in each target tissue, 
        based on ref_psi_target and measured delta_logit_psi in cat
        """
        if 'sample' not in self.df_mmsplice.columns:
            raise ValueError(
                '"sample" column is missing. Call add.samples() first')

        if type(cat_inference) == CatInference:
            cat_inference = [cat_inference]

        infer_rows = list()
        for cat in cat_inference:
            assert self.contains_chr == cat.contains_chr

            cols_shared = ['junction', 'gene_id', 'tissue', 'event_type']
            common_cat_idx = set(pd.concat([cat.common5, cat.common3])
                                 .set_index(['junctions', 'gene_id', 'tissue', 'event_type']).index)

            df_common = self.junction.reset_index().set_index(['junction', 'gene_id', 'tissue', 'event_type'])

            common_idx = set(df_common.index).intersection(common_cat_idx)
            df_common = df_common.loc[common_idx] \
                                 .set_index('sample', append=True)

            rows = df_common.index
            if progress:
                rows = tqdm(rows, total=df_common.shape[0])

            for junction, gene_id, tissue, event_type, sample in rows:
                if cat.contains(sample):
                    infer_rows.append(cat.infer(
                        junction, gene_id, tissue, sample, event_type))

        df = pd.DataFrame(infer_rows)
        assert df.shape[0] > 0
        df = df.set_index(['junction', 'gene_id', 'tissue', 'sample'])

        self.df_mmsplice_cat = self.junction.join(df)
        self.df_mmsplice_cat = self.df_mmsplice_cat[
            (~self.df_mmsplice_cat['tissue_cat'].isna())
            & (~self.df_mmsplice_cat['delta_psi_cat'].isna())
        ]

    def _get_maximum_effect(self, df, groupby, score):
        df = df.reset_index()
        if 'index' in df.columns:
            df = df.drop(columns='index')
        if len(set(groupby).difference(df.columns)) != 0:
            raise KeyError(" %s are not in columns" %
                           set(groupby).difference(df.columns))
        return get_abs_max_rows(df.set_index(groupby), groupby, score)

    @property
    def psi5(self):
        return SplicingOutlierResult(self.df_mmsplice[self.df_mmsplice['event_type'] == 'psi5'])

    @property
    def psi3(self):
        return SplicingOutlierResult(self.df_mmsplice[self.df_mmsplice['event_type'] == 'psi3'])

    @property
    def junction(self):  # NOTE: max aggregate over all variants
        groupby = ['junction', 'gene_id', 'event_type', 'tissue']
        if 'sample' in self.df_mmsplice:
            groupby.append('sample')
        if self._junction is None:
            self._junction = self._get_maximum_effect(
                self.df_mmsplice, groupby, score='delta_psi')
            self._junction = self._junction.reset_index('event_type')
        return self._junction

    @property
    def splice_site(self):  # NOTE: max aggregate over all variants
        groupby = ['splice_site', 'gene_id', 'event_type', 'tissue']
        if 'sample' in self.df_mmsplice:
            groupby.append('sample')
        if self._splice_site is None:
            self._splice_site = self._get_maximum_effect(
                self.df_mmsplice, groupby, score='delta_psi')
            self._splice_site = self._splice_site.reset_index('event_type')
        return self._splice_site

    @property
    def gene_mmsplice(self):  # NOTE: max aggregate over all variants
        groupby = ['gene_id', 'tissue']
        if 'sample' in self.df_mmsplice:
            groupby.append('sample')
        if self._gene_mmsplice is None:
            self._gene_mmsplice = self._get_maximum_effect(
                self.df_mmsplice, groupby, score='delta_psi')
        return self._gene_mmsplice

    @property
    # NOTE: max aggregate over all variants (and all CAT tissues)
    def gene_mmsplice_cat(self):
        # if wanted add tissue_cat in groupby
        groupby = ['gene_id', 'tissue', 'sample']
        if self._gene_mmsplice_cat is None:
            self._gene_mmsplice_cat = self._get_maximum_effect(
                self.df_mmsplice_cat, groupby, score='delta_psi_cat')
        return self._gene_mmsplice_cat

    @property
    def gene_spliceai(self):  # NOTE: max aggregate over all variants
        groupby = ['gene_id']
        if 'sample' in self.df_spliceai:
            groupby.append('sample')
        if self._gene_spliceai is None:
            self._gene_spliceai = self._get_maximum_effect(
                self.df_spliceai, groupby, score='delta_score')
        return self._gene_spliceai

    @property
    def variant_mmsplice(self):  # NOTE: max aggregate for variant on each gene
        groupby = ['variant', 'gene_id', 'tissue']
        if 'sample' in self.df_mmsplice:
            groupby.append('sample')
        if self._variant_mmsplice is None:
            self._variant_mmsplice = self._get_maximum_effect(
                self.df_mmsplice, groupby, score='delta_psi')
        return self._variant_mmsplice

    @property
    # NOTE: max aggregate for variant on each gene (over all CAT tissues)
    def variant_mmsplice_cat(self):
        # if wanted add tissue_cat in groupby
        groupby = ['variant', 'gene_id', 'tissue', 'sample']
        if self._variant_mmsplice_cat is None:
            self._variant_mmsplice_cat = self._get_maximum_effect(
                self.df_mmsplice_cat, groupby, score='delta_psi_cat')
        return self._variant_mmsplice_cat

    @property
    def variant_spliceai(self):  # NOTE: max aggregate for variant on each gene
        groupby = ['variant', 'gene_id']
        if 'sample' in self.df_spliceai:
            groupby.append('sample')
        if self._variant_spliceai is None:
            self._variant_spliceai = self._get_maximum_effect(
                self.df_spliceai, groupby, score='delta_score')
        return self._variant_spliceai

    @property
    def absplice_dna_input(self):
        if self._absplice_dna_input is None:
            groupby = ['variant', 'gene_id', 'tissue']
            if self.contains_samples:
                groupby.append('sample')
                
            # MMSplice (SpliceMap)
            cols_mmsplice = [
                'Chromosome', 'Start', 'End', 'Strand', 'junction', 'event_type', 'splice_site', 'gene_name',
                'delta_logit_psi', 'delta_psi', 'ref_psi', 'k', 'n', 'median_n',
                'novel_junction', 'weak_site_donor', 'weak_site_acceptor']
            if self.df_mmsplice is not None:
                df_mmsplice = self._get_maximum_effect(
                    self.df_mmsplice, groupby, score='delta_psi')
            else:
                df_mmsplice = pd.DataFrame(columns=[*cols_mmsplice, *groupby]).set_index(groupby)
            
            # SpliceAI
            cols_spliceai = ['delta_score', 'gene_name']
            if self.df_spliceai is not None:
                df_spliceai = self._add_tissue_info_to_spliceai()
                df_spliceai = self._get_maximum_effect(
                    df_spliceai, groupby, score='delta_score')
            else:
                df_spliceai = pd.DataFrame(columns=[*cols_spliceai, *groupby]).set_index(groupby)

            # Join MMSplice & SpliceAI
            self._absplice_dna_input = df_mmsplice[cols_mmsplice].join(
                df_spliceai[cols_spliceai], how='outer', rsuffix='_spliceai')

            self._absplice_dna_input = self._absplice_dna_input.reset_index()
            self._absplice_dna_input = self._absplice_dna_input.set_index(['gene_id', 'tissue'])\
                .join(self.gene_tpm.set_index(['gene_id', 'tissue'])[['gene_tpm']])\
                .reset_index().set_index(groupby)
        return self._absplice_dna_input

    @property
    def absplice_rna_input(self):
        if self._absplice_rna_input is None:
            groupby = ['variant', 'gene_id', 'tissue', 'sample']
            if not pd.Series(groupby).isin(self.absplice_dna_input.index.names).all():
                self._absplice_dna_input = self.absplice_dna_input.set_index(
                    groupby)
            df_mmsplice_cat = self._get_maximum_effect(
                self.df_mmsplice_cat, groupby, score='delta_psi_cat')
            cols_mmsplice_cat = [
                'junction', 'delta_psi', 'ref_psi', 'median_n',
                *[col for col in df_mmsplice_cat.columns if 'cat' in col]]
            self._absplice_rna_input = self.absplice_dna_input.join(
                df_mmsplice_cat[cols_mmsplice_cat], how='outer', rsuffix='_from_cat_infer')
        return self._absplice_rna_input

    def _predict_absplice(self, df, absplice_score, pickle_file, features, abs_features, median_n_cutoff, tpm_cutoff):
        model = pickle.load(open(pickle_file, 'rb'))
        df['splice_site_is_expressed'] = (
            df['median_n'] > median_n_cutoff).astype(int)
        df['gene_is_expressed'] = (df['gene_tpm'] > tpm_cutoff).astype(int)
        df = df[features].fillna(0)
        if abs_features:
            df = np.abs(df)
        df[absplice_score] = model.predict_proba(df)[:, 1]
        return df

    def predict_absplice_dna(self, pickle_file=None, features=None, abs_features=False, median_n_cutoff=0, tpm_cutoff=1):
        if features is None:
            features = [
                'delta_logit_psi',
                'delta_psi', 
                'delta_score',
                'splice_site_is_expressed']
        if pickle_file is None:
            pickle_file = ABSPLICE_DNA

        self._absplice_dna = self._predict_absplice(
            df=self.absplice_dna_input,
            absplice_score='AbSplice_DNA',
            pickle_file=pickle_file,
            features=features,
            abs_features=abs_features,
            median_n_cutoff=median_n_cutoff,
            tpm_cutoff=tpm_cutoff)
        return self._absplice_dna

    def predict_absplice_rna(self, pickle_file=None, features=None, abs_features=False, median_n_cutoff=0, tpm_cutoff=1):
        if features is None:
            features = [
                'delta_logit_psi',
                'delta_psi',
                'delta_psi_cat',
                'delta_score',
                'splice_site_is_expressed']
        if pickle_file is None:
            pickle_file = ABSPLICE_RNA

        self._absplice_rna = self._predict_absplice(
            df=self.absplice_rna_input,
            absplice_score='AbSplice_RNA',
            pickle_file=pickle_file,
            features=features,
            abs_features=abs_features,
            median_n_cutoff=median_n_cutoff,
            tpm_cutoff=tpm_cutoff)
        return self._absplice_rna

    @property
    def gene_absplice_dna(self):  # NOTE: max aggregate over all variants
        groupby = ['gene_id', 'tissue']
        if 'sample' in self._absplice_dna.reset_index():
            groupby.append('sample')
        if self._gene_absplice_dna is None:
            self._gene_absplice_dna = self._get_maximum_effect(
                self._absplice_dna, groupby, score='AbSplice_DNA')
        return self._gene_absplice_dna

    @property
    def gene_absplice_rna(self):  # NOTE: max aggregate over all variants
        groupby = ['gene_id', 'tissue', 'sample']
        if self._gene_absplice_rna is None:
            self._gene_absplice_rna = self._get_maximum_effect(
                self._absplice_rna, groupby, score='AbSplice_RNA')
        return self._gene_absplice_rna

    @property
    # NOTE: max aggregate for variant on each gene
    def variant_absplice_dna(self):
        groupby = ['variant', 'gene_id', 'tissue']
        if 'sample' in self._absplice_dna.reset_index():
            groupby.append('sample')
        if self._variant_absplice_dna is None:
            self._variant_absplice_dna = self._get_maximum_effect(
                self._absplice_dna, groupby, score='AbSplice_DNA')
        return self._variant_absplice_dna

    @property
    # NOTE: max aggregate for variant on each gene
    def variant_absplice_rna(self):
        groupby = ['variant', 'gene_id', 'tissue', 'sample']
        if self._variant_absplice_rna is None:
            self._variant_absplice_rna = self._get_maximum_effect(
                self._absplice_rna, groupby, score='AbSplice_RNA')
        return self._variant_absplice_rna

    @staticmethod
    def _add_maf(df, population, default=-1):
        df['maf'] = df['variant'].map(
            lambda x: population.get(x, default))
        return df

    @staticmethod
    def _filter_private(df, max_num_sample=2):
        df['samples'] = df.groupby(
            'variant')['sample'].apply(lambda x: ';'.join(x))
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
            if df_mmsplice is not None:
                df_mmsplice = self._filter_private(df_mmsplice, max_num_sample)
            if df_spliceai is not None:
                df_spliceai = self._filter_private(df_spliceai, max_num_sample)

        if population:
            if df_mmsplice is not None:
                df_mmsplice = self._add_filter_maf(
                    df_mmsplice, population, maf_cutoff, default)
            if df_spliceai is not None:
                df_spliceai = self._add_filter_maf(
                    df_spliceai, population, maf_cutoff, default)

        return SplicingOutlierResult(
            df_mmsplice=df_mmsplice,
            df_spliceai=df_spliceai
        )