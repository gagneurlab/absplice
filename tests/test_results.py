import pytest
import pandas as pd
import numpy as np
from kipoiseq.extractors.vcf import MultiSampleVCF
# from kipoiseq.extractors.vcf_query import to_sample_csv
from absplice import SpliceOutlier, SpliceOutlierDataloader, CatInference, SplicingOutlierResult
from absplice.ensemble import train_model_ebm
from absplice.utils import inject_new_row
from absplice.result import GENE_MAP, GENE_TPM
from conftest import df_mmsplice_cat, multi_vcf_file, \
    mmsplice_path, spliceai_path, mmsplice_cat_path, var_samples_path, \
        fasta_file, ref_table5_kn_testis, ref_table5_kn_lung, ref_table3_kn_testis, ref_table3_kn_lung, \
            spliceai_vcf_path2, cadd_splice_path, outliers_cat_path
        # absplice_precomputed_path
    
def test_splicing_outlier_result__init__mmsplice_only(df_mmsplice):
    # initialize with pd.DataFrame
    sor = SplicingOutlierResult(
        df_mmsplice = df_mmsplice
    )
    assert sor.df_mmsplice is not None
    assert 'delta_psi' in sor.df_mmsplice.columns
    
    # initialitze with path
    sor = SplicingOutlierResult(
        df_mmsplice = mmsplice_path
    )
    assert sor.df_mmsplice is not None
    assert 'delta_psi' in sor.df_mmsplice.columns
    
def test_splicing_outlier_result__init__spliceai_only(df_spliceai):
    # initialize with pd.DataFrame
    sor = SplicingOutlierResult(
        df_spliceai = df_spliceai
    )
    assert sor.df_spliceai is not None
    assert 'delta_score' in sor.df_spliceai.columns
    del sor
    
    # initialitze with path
    sor = SplicingOutlierResult(
        df_spliceai = spliceai_path
    )
    assert sor.df_spliceai is not None
    assert 'delta_score' in sor.df_spliceai.columns
    del sor
    
    # initialitze with vcf
    sor = SplicingOutlierResult(
        df_spliceai = spliceai_vcf_path2
    )
    assert sor.df_spliceai is not None
    assert 'delta_score' in sor.df_spliceai.columns
    
    
def test_splicing_outlier_result__init__cadd_splice_only():
    # # initialize with pd.DataFrame
    # sor = SplicingOutlierResult(
    #     df_cadd_splice = df_cadd_splice
    # )
    # assert sor.df_cadd_splice is not None
    # assert 'PHRED' in sor.df_cadd_splice.columns
    
    # initialitze with path
    sor = SplicingOutlierResult(
        df_cadd_splice = cadd_splice_path
    )
    assert sor.df_cadd_splice is not None
    assert 'PHRED' in sor.df_cadd_splice.columns
    
    
def test_splicing_outlier_result__init__mmsplice_cat_only(df_mmsplice_cat):
    # initialize with pd.DataFrame
    sor = SplicingOutlierResult(
        df_mmsplice_cat = df_mmsplice_cat
    )
    assert sor.df_mmsplice_cat is not None
    assert 'delta_psi_cat' in sor.df_mmsplice_cat.columns
    
    # initialitze with path
    sor = SplicingOutlierResult(
        df_mmsplice_cat = mmsplice_cat_path
    )
    assert sor.df_mmsplice_cat is not None
    assert 'delta_psi_cat' in sor.df_mmsplice_cat.columns
    
def test_splicing_outlier_result__init__absplice_dna_input():
    sor_absplice_dna = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path, 
        gene_tpm=GENE_TPM,
        gene_map=GENE_MAP
    )
    df_absplice_dna_input = sor_absplice_dna.absplice_dna_input
    
    sor = SplicingOutlierResult(
        df_absplice_dna_input=df_absplice_dna_input
    )
    assert sor.absplice_dna_input.shape[0] > 0
    
    
# def test_splicing_outlier_result__init__absplice_dna_input_precomputed():
#     sor = SplicingOutlierResult(
#         df_absplice_dna_input = absplice_precomputed_path
#     )
#     assert sor.absplice_dna_input.shape[0] > 0
    
    
def test_splicing_outlier_result__init__absplice_dna_input_CADD():
    sor_absplice_dna = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path, 
        df_cadd_splice=cadd_splice_path,
        gene_tpm=GENE_TPM,
        gene_map=GENE_MAP
    )
    df_absplice_dna_input = sor_absplice_dna.absplice_dna_input
    
    sor = SplicingOutlierResult(
        df_absplice_dna_input=df_absplice_dna_input
    )
    assert sor.absplice_dna_input.shape[0] > 0
    assert 'PHRED' in sor.absplice_dna_input.columns
    
    
def test_splicing_outlier_result__init__absplice_dna_input_mmsplice_None():
    sor_absplice_dna = SplicingOutlierResult(
        df_mmsplice=None, 
        df_spliceai=spliceai_path, 
        gene_tpm=GENE_TPM,
        gene_map=GENE_MAP
    )
    df_absplice_dna_input = sor_absplice_dna.absplice_dna_input
    
    sor = SplicingOutlierResult(
        df_absplice_dna_input=df_absplice_dna_input
    )
    assert sor.absplice_dna_input.shape[0] > 0
    
    
def test_splicing_outlier_result__init__absplice_dna_input_spliceai_None():
    sor_absplice_dna = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=None, 
        gene_tpm=GENE_TPM,
        gene_map=GENE_MAP
    )
    df_absplice_dna_input = sor_absplice_dna.absplice_dna_input
    
    sor = SplicingOutlierResult(
        df_absplice_dna_input=df_absplice_dna_input
    )
    assert sor.absplice_dna_input.shape[0] > 0
    
    
def test_splicing_outlier_result__init__absplice_rna_input():
    # Initialize with mmsplice_cat
    sor_absplice_rna = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path, 
        df_mmsplice_cat=mmsplice_cat_path, 
        df_outliers_cat=outliers_cat_path,
        gene_tpm=GENE_TPM,
        gene_map=GENE_MAP,
        df_var_samples=var_samples_path
    )
    df_absplice_rna_input = sor_absplice_rna.absplice_rna_input
    sor = SplicingOutlierResult(
        df_absplice_rna_input=df_absplice_rna_input
    )
    assert sor.absplice_rna_input.shape[0] > 0
    
    
# def test_splicing_outlier_result__init__absplice_rna_input_absplice_dna_precomputed():
#     # Initialize with mmsplice_cat
#     sor_absplice_rna = SplicingOutlierResult(
#         df_absplice_dna_input = absplice_precomputed_path,
#         df_mmsplice_cat=mmsplice_cat_path, 
#         df_var_samples=var_samples_path
#     )
#     df_absplice_rna_input = sor_absplice_rna.absplice_rna_input
#     sor = SplicingOutlierResult(
#         df_absplice_rna_input=df_absplice_rna_input
#     )
#     assert sor.absplice_rna_input.shape[0] > 0
     
     
     
# def test_splicing_outlier_result_add_spliceai(outlier_results, df_spliceai, gene_map):
#     assert outlier_results.df_spliceai is None
#     outlier_results.add_spliceai(df_spliceai, gene_map)
#     assert outlier_results.df_spliceai is not None
    
# def test_splicing_outlier_result_add_tissue_info_to_spliceai(outlier_results_complete):
#     results = outlier_results_complete
#     results._add_tissue_info_to_spliceai()
#     assert 'tissue' in results._df_spliceai_tissue.columns
    
def test_splicing_outlier_result_add_samples():
    sor = SplicingOutlierResult(
        df_mmsplice = mmsplice_path,
        df_spliceai = spliceai_path,
        gene_map = GENE_MAP
    )
    assert 'sample' not in sor.df_mmsplice.columns
    assert 'sample' not in sor.df_spliceai.columns
    
    sor.add_samples(var_samples_path)
    assert 'sample' in sor.df_mmsplice.columns
    assert 'sample' in sor.df_spliceai.columns
    
def test_splicing_outlier_result_init_add_samples():
    sor = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path, 
        df_var_samples = var_samples_path,
    )
    
    sor2 = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path, 
    )
    sor2.add_samples(var_samples_path)
    
    assert sor.df_mmsplice.equals(sor2.df_mmsplice)
    assert sor.df_spliceai.equals(sor2.df_spliceai)
    
def test_splicing_outlier_result_infer_cat(outlier_results_multi, cat_dl, df_var_samples):
    
    results = outlier_results_multi
    results.add_samples(df_var_samples)
    results.infer_cat(cat_dl)

    assert sorted(results.df_mmsplice_cat.reset_index().columns.tolist()) == sorted([
        'variant', 'tissue', 'junction', 'event_type',
        'splice_site', 'ref_psi', 'median_n', 
        'gene_id', 'gene_name',
        'delta_logit_psi', 'delta_psi',
        'sample', 'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
        'delta_logit_psi_cat', 'delta_psi_cat'])

    assert results.df_mmsplice_cat.loc[(
        '17:41201211-41203079:-', 'ENSG00000012048', 'Testis', 'NA00002',)] is not None
    
def test_splicing_outlier_result_psi(outlier_results):
    results = outlier_results.psi5
    assert all(results.psi5.df_mmsplice['event_type'] == 'psi5')

    results = outlier_results.psi3
    assert all(results.psi3.df_mmsplice['event_type'] == 'psi3')

def test_splicing_outlier_result_splice_site(outlier_results):
    assert sorted(outlier_results.splice_site.index) \
        == sorted(set(outlier_results.df_mmsplice.set_index(['splice_site', 'gene_id', 'tissue']).index))
    
def test_splicing_outlier_result_gene_mmsplice(outlier_results, outlier_results_multi):
    assert sorted(set(outlier_results.gene_mmsplice.index)) \
        == sorted(set(outlier_results.df_mmsplice.set_index(['gene_id', 'tissue']).index))

    assert sorted(set(outlier_results_multi.gene_mmsplice.index)) \
        == sorted(set(outlier_results_multi.df_mmsplice.set_index(['gene_id', 'tissue', 'sample']).index))
       
def test_splicing_outlier_result_gene_mmsplice_max_effect():
    sor = SplicingOutlierResult(
        df_mmsplice = mmsplice_path
    )
    
    assert sor.df_mmsplice.groupby('gene_id')['delta_psi'].apply(
        lambda x: np.abs(x).max())['ENSG00000012048'] == 0.1023274783099152
    assert np.abs(sor.gene_mmsplice['delta_psi'].values[0]) == 0.1023274783099152
    
    # inject new variant with higher effect
    sor._gene_mmsplice = None  
    sor.df_mmsplice = inject_new_row(sor.df_mmsplice, new_row_dict = {
        'variant': '1:1:T>C',
        'delta_psi': 1.0})
    assert np.abs(sor.gene_mmsplice['delta_psi'].values[0]) == 1.0
    
    # inject new gene_id
    sor._gene_mmsplice = None
    sor.df_mmsplice = inject_new_row(sor.df_mmsplice, new_row_dict = {
        'gene_id': 'test_gene',
        'variant': '1:1:T>C',
        'delta_psi': 1.0})
    assert 'test_gene' in sor.gene_mmsplice.index.get_level_values('gene_id').values
     
def test_splicing_outlier_result_gene_mmsplice_cat(outlier_dl_multi, cat_dl, outlier_model, df_var_samples):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.add_samples(df_var_samples)
    results.infer_cat(cat_dl)
    
    assert results.df_mmsplice_cat is not None
    assert 'tissue_cat' in results.df_mmsplice_cat.columns
    
    # all missing junctions do not have variant in vicinity (no mmsplice predictions either)
    df = results.df_mmsplice_cat[results.df_mmsplice_cat['tissue_cat'] == cat_dl[0].ct.name]
    common_junctions = cat_dl[0].common_junctions3[0].union(cat_dl[0].common_junctions5[0])
    assert len(set(results.df_mmsplice['junction'])\
        .intersection(common_junctions.difference(df.index.get_level_values('junction')))) == 0
    
def test_splicing_outlier_result_variant_mmsplice_max_effect():
    sor = SplicingOutlierResult(
        df_mmsplice = mmsplice_path
    )

    assert np.abs(sor.variant_mmsplice['delta_psi'].values[0]) == 0.1023274783099152
    
    # inject new variant with higher effect
    sor._variant_mmsplice = None  
    sor.df_mmsplice = inject_new_row(sor.df_mmsplice, new_row_dict = {
        'variant': '1:1:T>C',
        'delta_psi': 1.0})
    assert sor.variant_mmsplice[['delta_psi']].values.max() == 1.0
    
    # inject new variant on different gene_id
    sor._variant_mmsplice = None
    sor.df_mmsplice = inject_new_row(sor.df_mmsplice, new_row_dict = {
        'gene_id': 'test_gene',
        'variant': '1:1:T>C',
        'delta_psi': 1.0})
    assert 'test_gene' in sor.variant_mmsplice.index.get_level_values('gene_id').values
    
    
def test_splicing_outlier_result_predict_absplice_dna():
    sor = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path,
    )
    sor.predict_absplice_dna()
    assert 'AbSplice_DNA' in sor._absplice_dna.columns
    
    
def test_splicing_outlier_result_predict_absplice_dna_CADD():
    sor = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path,
        df_cadd_splice=cadd_splice_path
    )
    sor.predict_absplice_dna(cadd_splice=True)
    assert 'AbSplice_DNA' in sor._absplice_dna.columns
    assert 'PHRED' in sor._absplice_dna.columns
    
    
def test_splicing_outlier_result_predict_absplice_rna():
    sor = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path,
        df_mmsplice_cat=mmsplice_cat_path,
        df_outliers_cat=outliers_cat_path,
        df_var_samples=var_samples_path
    )
    sor.predict_absplice_rna()
    assert 'AbSplice_RNA' in sor._absplice_rna.columns


def test_result_absplice_dna_sorting():
    df_mmsplice = pd.DataFrame({
        'variant': ['v1', 'v1'],
        'gene_id': ['g1', 'g1'],
        'gene_name': ['gene1', 'gene1'],
        'event_type': ['psi3', 'psi3'],
        'splice_site': ['chr:1:-', 'chr:1:-'],
        'tissue': ['t1', 't1'],
        'junction': ['j1', 'j2'],
        'delta_psi': [0.1, 0.101],
        'delta_logit_psi': [0.1, 0.1],
        'ref_psi': [0.4, 0.4],
        'median_n': [100, 2],
    })

    df_spliceai = pd.DataFrame({
        'variant': ['v1'],
        'gene_name': ['gene1'],
        'delta_score': [0.12],
        'acceptor_gain': [0.1],
        'acceptor_loss': [0.12],
        'donor_gain': [0.1],
        'donor_loss': [0.1],
        'acceptor_gain_position': [2],
        'acceptor_loss_position': [2],
        'donor_gain_position': [2],
        'donor_loss_position': [2]
    })

    splicing_result = SplicingOutlierResult(
            df_mmsplice=df_mmsplice, 
            df_spliceai=df_spliceai,
            gene_map = pd.DataFrame({'gene_id': ['g1'], 'gene_name': ['gene1']})
        )
    
    splicing_result.predict_absplice_dna()
    assert splicing_result._absplice_dna.shape == (1,18)
    assert splicing_result._absplice_dna['junction'].values == ['j1']

    df_mmsplice['sample'] = ['s1', 's1']
    df_spliceai['sample'] = ['s1']

    splicing_result = SplicingOutlierResult(
        df_mmsplice=df_mmsplice, 
        df_spliceai=df_spliceai,
        gene_map = pd.DataFrame({'gene_id': ['g1'], 'gene_name': ['gene1']})
    )
    
    splicing_result.predict_absplice_dna()
    assert splicing_result._absplice_dna.shape == (1,18)
    assert len(splicing_result._absplice_dna.index[0]) == 4
    assert splicing_result._absplice_dna['junction'].values == ['j1']

def test_absplice_rna_input_join():
    d1 = pd.DataFrame({'variant': ['v1', 'v1'], 'score': [0.2, 0]})
    d2 = pd.DataFrame({'variant': ['v1', 'v1'], 'score': [0.3, 0.2], 'tissue_cat': ['tc1', 'tc2']})


    df = d1.set_index('variant').join(d2.set_index('variant'), how='outer', rsuffix='_cat')

    assert df.shape == (4,3)
    assert True not in df['tissue_cat'].isna().values

def test_absplice_rna_output_cat_difference():
    df_rna_input_test = pd.DataFrame({
        'variant': ['v1', 'v1', 'v1'],
        'gene_id': ['g1', 'g1', 'g1'],
        'tissue': ['t1', 't1', 't1'],
        'sample': ['s1', 's1', 's1'],
        'junction': ['j1', 'j1', 'j2'],
        'event_type': ['psi3', 'psi3', 'psi5'],
        'splice_site': ['chr:1:-', 'chr:1:-', 'chr:1:-'],
        'delta_score': [0.5, 0.5, 0.2],
        'acceptor_gain': [0.5, 0.5, 0.1],
        'acceptor_loss': [0.1, 0.1, 0.2],
        'donor_gain': [0,0,0],
        'donor_loss': [0,0,0],
        'acceptor_gain_position': [2,2,3],
        'acceptor_loss_position': [2,2,3],
        'donor_gain_position': [2,2,3],
        'donor_loss_position': [2,2,3],
        'delta_logit_psi': [0.7, 0.7, 0.4],
        'delta_psi': [0.4, 0.4, 0.1], 
        'ref_psi': [0.2, 0.2, 0.2], 
        'median_n': [30, 30, 5],
        'tissue_cat': ['tc1', 'tc2', 'tc1'], 
        'k_cat': [10, 10, 10],
        'n_cat': [10, 10, 10], 
        'median_n_cat': [60, 100, 12], 
        'psi_cat': [0.7, 0.2, 0.7], 
        'ref_psi_cat': [0.1, 0.2, 0.4],
        'delta_logit_psi_cat': [0.9, 0, 0.6], 
        'delta_psi_cat': [0.6, 0, 0.3], 
        'pValueGene_g_minus_log10': [0,0,0]
    })

    splicing_result = SplicingOutlierResult(
        df_absplice_rna_input=df_rna_input_test,
        gene_map = pd.DataFrame({'gene_id': ['g1'], 'gene_name': ['gene1']})
    )

    splicing_result.predict_absplice_rna()

    assert splicing_result._absplice_rna.shape == (1,19)
    assert splicing_result._absplice_rna['delta_psi_cat'].values[0] == 0.6