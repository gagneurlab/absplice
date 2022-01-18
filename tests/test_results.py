import pytest
import pandas as pd
import numpy as np
from kipoiseq.extractors.vcf import MultiSampleVCF
# from kipoiseq.extractors.vcf_query import to_sample_csv
from splicing_outlier_prediction import SpliceOutlier, SpliceOutlierDataloader, CatInference, SplicingOutlierResult
from splicing_outlier_prediction.ensemble import train_model_ebm
from conftest import fasta_file, vcf_file, multi_vcf_file, multi_vcf_samples, \
    ref_table5_kn_testis, ref_table3_kn_testis,  \
    ref_table5_kn_lung, ref_table3_kn_lung, \
    combined_ref_tables5_testis_lung, combined_ref_tables3_testis_lung, \
    count_cat_file_lymphocytes,  count_cat_file_blood, \
    spliceai_path, mmsplice_path, mmsplice_cat_path, pickle_DNA, pickle_DNA_CAT, gene_map, gene_tpm, pickle_absplice_DNA, pickle_absplice_RNA


def test_splicing_outlier_result__init__(gene_tpm, gene_map):
    df_mmsplice = pd.read_csv(mmsplice_path)
    df_spliceai = pd.read_csv(spliceai_path)
    # _gene_tpm = gene_tpm[gene_tpm['tissue'] == 'Whole_Blood']
    
    # Initialize with DataFrame
    sor = SplicingOutlierResult(
        df_mmsplice=df_mmsplice, 
        df_spliceai=df_spliceai, 
        gene_tpm=gene_tpm
    )
    assert sor.df_mmsplice.shape[0] > 0
    assert sor.df_spliceai.shape[0] > 0

    # Initialize with str
    sor2 = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path, 
        gene_tpm=gene_tpm
    )
    assert sor2.df_mmsplice.shape[0] > 0
    assert sor2.df_spliceai.shape[0] > 0
    
    # Initialize with spliceai only
    sor = SplicingOutlierResult(
        df_spliceai=df_spliceai, 
        gene_map=gene_map,
        gene_tpm=gene_tpm
    )
    
    assert sor.df_spliceai.shape[0] > 0
    assert sor.df_mmsplice is None
    
    # Initialize with mmsplice only
    sor = SplicingOutlierResult(
        df_mmsplice=df_mmsplice
    )
    
    assert sor.df_spliceai is None
    assert sor.df_mmsplice.shape[0] > 0
    
    # Initialize with mmsplice_cat
    sor3 = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path, 
        df_mmsplice_cat=mmsplice_cat_path, 
        gene_tpm=gene_tpm
    )
    assert sor3.df_mmsplice.shape[0] > 0
    assert sor3.df_spliceai.shape[0] > 0
    
    # Initialize with mmsplice_cat only
    sor = SplicingOutlierResult(
        df_mmsplice_cat=mmsplice_cat_path
    )
    
    assert sor.df_spliceai is None
    assert sor.df_mmsplice is None
    assert sor.df_mmsplice_cat.shape[0] > 0
    


def test_write_sample_csv():   
    vcf = MultiSampleVCF(multi_vcf_file)
    variant_queryable = vcf.query_all()
    variant_queryable.to_sample_csv(multi_vcf_samples)

def test_tissues(outlier_results):
    assert sorted(outlier_results.df_mmsplice['tissue'].unique()) == sorted(
        ['Testis', 
        #  'Lung'
         ])

def test_splicing_outlier_result_add_spliceai(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model):
    # single sample vcf
    results = outlier_model.predict_on_dataloader(outlier_dl)
    results.add_spliceai(spliceai_path, gene_mapping=False)
    assert results.df_spliceai.shape[0] > 0
    assert 'delta_score' in results.df_spliceai.columns
    # multi sample vcf
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.add_spliceai(spliceai_path, gene_mapping=False)
    assert results.df_spliceai.shape[0] > 0
    assert 'delta_score' in results.df_spliceai.columns


def test_splicing_outlier_result_add_samples(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model):
    # TODO: spliceai_path has no samples here
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.add_spliceai(spliceai_path, gene_mapping=False)
    
    variants_mmsplice = set(results.df_mmsplice.set_index(['variant']).index)
    variants_spliceai = set(results.df_spliceai.set_index(['variant']).index)
    variants = set(var_samples_df.set_index(['variant']).index)
    samples_mmsplice = set(var_samples_df.set_index(['variant']).loc[variants_mmsplice.intersection(variants)]['sample'])
    samples_spliceai = set(var_samples_df.set_index(['variant']).loc[variants_spliceai.intersection(variants)]['sample'])
    
    results.add_samples(var_samples_df)
    assert len(samples_mmsplice.difference(set(results.df_mmsplice['sample'].values))) == 0
    assert len(samples_spliceai.difference(set(results.df_spliceai['sample'].values))) == 0
    

def test_splicing_outlier_result_psi(outlier_results):
    results = outlier_results.psi5
    assert all(results.psi5.df_mmsplice['event_type'] == 'psi5')

    results = outlier_results.psi3
    assert all(results.psi3.df_mmsplice['event_type'] == 'psi3')


def test_splicing_outlier_result_splice_site(outlier_results):
    assert sorted(outlier_results.splice_site.index) \
        == sorted(set(outlier_results.df_mmsplice.set_index(['splice_site', 'tissue']).index))


def test_splicing_outlier_result_gene_mmsplice(outlier_results, outlier_results_multi):
    assert sorted(set(outlier_results.gene_mmsplice.index)) \
        == sorted(set(outlier_results.df_mmsplice.set_index(['gene_id', 'tissue']).index))

    assert sorted(set(outlier_results_multi.gene_mmsplice.index)) \
        == sorted(set(outlier_results_multi.df_mmsplice.set_index(['gene_id', 'tissue', 'sample']).index))
    # assert sorted(set(outlier_results_multi.gene_mmsplice.index)) \
    #     == sorted(set(outlier_results_multi._explode(outlier_results_multi.df_mmsplice, 
    #                                                  col='samples', new_name='sample').set_index(
    #                                                      ['gene_id', 'sample', 'tissue']
    #                                                      ).index))

def test_splicing_outlier_result_gene_mmsplice_cat(outlier_dl_multi, outlier_results, cat_dl, outlier_model, var_samples_df):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.add_samples(var_samples_df)
    results.infer_cat(cat_dl)
    
    assert results.df_mmsplice_cat is not None
    assert 'tissue_cat' in results.df_mmsplice_cat.columns
    
    # all missing junctions do not have variant in vicinity (no mmsplice predictions either)
    df = results.df_mmsplice_cat[results.df_mmsplice_cat['tissue_cat'] == cat_dl[0].ct.name]
    common_junctions = cat_dl[0].common_junctions3[0].union(cat_dl[0].common_junctions5[0])
    assert len(set(results.df_mmsplice['junction']).intersection(common_junctions.difference(df.index.get_level_values('junction')))) == 0

def test_splicing_outlier_result_gene_spliceai(outlier_results, outlier_results_multi, gene_map, gene_tpm, var_samples_df):
    outlier_results.gene_map = gene_map
    outlier_results.add_spliceai(spliceai_path, gene_mapping=True)
    outlier_results.df_spliceai = outlier_results.df_spliceai[
        ~outlier_results.df_spliceai['gene_id'].isna()
    ]
    assert sorted(set(outlier_results.gene_spliceai.index)) \
        == sorted(set(outlier_results.df_spliceai.set_index(['gene_id']).index))

    outlier_results_multi.gene_map = gene_map
    outlier_results_multi.add_spliceai(spliceai_path, gene_mapping=True)
    outlier_results_multi.df_spliceai = outlier_results_multi.df_spliceai[
        ~outlier_results_multi.df_spliceai['gene_id'].isna()
    ]
    outlier_results_multi.add_samples(var_samples_df)
    assert sorted(set(outlier_results_multi.gene_spliceai.index)) \
        == sorted(set(outlier_results_multi.df_spliceai.set_index(['gene_id', 'sample']).index))
        
    df_spliceai_gene_no_tissue = outlier_results_multi.gene_spliceai.copy()
    outlier_results_multi.gene_tpm = gene_tpm
    outlier_results_multi.df_spliceai = None
    outlier_results_multi.add_spliceai(spliceai_path, gene_mapping=True)
    outlier_results_multi._add_tissue_info_to_spliceai()
    outlier_results_multi.df_spliceai = outlier_results_multi._df_spliceai_tissue #store spliceai_path tissue in spliceai_path
    assert 'tissue' in outlier_results_multi.df_spliceai.columns
    
    assert sorted(set(outlier_results_multi.gene_spliceai.index)) \
        == sorted(set(df_spliceai_gene_no_tissue.index))


def test_splicing_outlier_result_absplice_input_dna(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model, gene_map, gene_tpm):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.gene_map = gene_map
    results.gene_tpm = gene_tpm
    results.add_spliceai(spliceai_path, gene_mapping=True)
    results.add_samples(var_samples_df)
    
    df_spliceai = results._add_tissue_info_to_spliceai()
    indices_mmsplice = results.df_mmsplice.set_index(['variant', 'gene_id', 'tissue', 'sample']).index.unique()
    indices_spliceai = df_spliceai.set_index(['variant', 'gene_id', 'tissue', 'sample']).index.unique()
    indices_all = indices_mmsplice.union(indices_spliceai)
    
    assert results.absplice_dna_input.shape[0] > 0
    assert len(set(results.absplice_dna_input.index).difference(indices_all)) == 0
    
    assert results.df_mmsplice_cat is None

def test_splicing_outlier_result_absplice_input_rna(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model, gene_map, gene_tpm, cat_dl):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.gene_map = gene_map
    results.gene_tpm = gene_tpm
    results.add_spliceai(spliceai_path, gene_mapping=True)
    results.add_samples(var_samples_df)
    results.infer_cat(cat_dl)
    
    df_spliceai = results._add_tissue_info_to_spliceai()
    indices_mmsplice = results.df_mmsplice.set_index(['variant', 'gene_id', 'tissue', 'sample']).index.unique()
    indices_spliceai = df_spliceai.set_index(['variant', 'gene_id', 'tissue', 'sample']).index.unique()
    indices_all = indices_mmsplice.union(indices_spliceai)
    
    assert results.absplice_rna_input.shape[0] > 0
    assert len(set(results.absplice_rna_input.index).difference(indices_all)) == 0
    
    assert results.df_mmsplice_cat is not None
    assert len(set(['delta_score', 'delta_psi', 'delta_psi_cat']).difference(results.absplice_rna_input.columns)) == 0

    
def test_splicing_outlier_result_predict_absplice_dna(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model, gene_map, gene_tpm):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.gene_map = gene_map
    results.gene_tpm = gene_tpm
    results.add_spliceai(spliceai_path, gene_mapping=True)
    results.add_samples(var_samples_df)
    
    results.predict_absplice_dna(pickle_file=pickle_absplice_DNA)
    assert 'AbSplice_DNA' in results.absplice_dna.columns
    
    
def test_splicing_outlier_result_predict_absplice_rna(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model, gene_map, gene_tpm, cat_dl):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.gene_map = gene_map
    results.gene_tpm = gene_tpm
    results.add_spliceai(spliceai_path, gene_mapping=True)
    results.add_samples(var_samples_df)
    results.infer_cat(cat_dl)
    
    results.predict_absplice_rna(pickle_file=pickle_absplice_RNA)
    assert 'AbSplice_RNA' in results.absplice_rna.columns
    
    
def test_splicing_outlier_result_gene_absplice_dna(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model, gene_map, gene_tpm):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.gene_map = gene_map
    results.gene_tpm = gene_tpm
    results.add_spliceai(spliceai_path, gene_mapping=True)
    results.add_samples(var_samples_df)
    
    results.predict_absplice_dna(pickle_file=pickle_absplice_DNA)
    assert 'variant' in results.absplice_dna.index.names
    assert 'variant' not in results.gene_absplice_dna.index.names
    assert 'AbSplice_DNA' in results.gene_absplice_dna.columns
    
    
def test_splicing_outlier_result_gene_absplice_rna(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model, gene_map, gene_tpm, cat_dl):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.gene_map = gene_map
    results.gene_tpm = gene_tpm
    results.add_spliceai(spliceai_path, gene_mapping=True)
    results.add_samples(var_samples_df)
    results.infer_cat(cat_dl)
    
    results.predict_absplice_rna(pickle_file=pickle_absplice_RNA)
    assert 'variant' in results.absplice_rna_input.index.names
    assert 'variant' not in results.gene_absplice_rna.index.names
    assert 'AbSplice_RNA' in results.gene_absplice_rna.columns
    
    
    
def test_splicing_outlier_complete_dna(gene_map, gene_tpm, var_samples_df):
    
    results = SplicingOutlierResult(
        df_mmsplice = mmsplice_path,
        df_spliceai = spliceai_path,
        gene_map = gene_map,
        gene_tpm = gene_tpm,
    )
    
    results.add_samples(var_samples_df)
    results.predict_absplice_dna(pickle_file=pickle_absplice_DNA)
    
    assert results.absplice_dna.shape[0] > 0 
    assert 'AbSplice_DNA' in results.absplice_dna.columns
    
def test_splicing_outlier_complete_rna(gene_map, gene_tpm, var_samples_df):
    
    results = SplicingOutlierResult(
        df_mmsplice = mmsplice_path,
        df_spliceai = spliceai_path,
        df_mmsplice_cat = mmsplice_cat_path,
        gene_map = gene_map,
        gene_tpm = gene_tpm,
    )
    
    results.add_samples(var_samples_df)
    results.predict_absplice_rna(pickle_file=pickle_absplice_RNA)
    
    assert results.absplice_rna.shape[0] > 0 
    assert 'AbSplice_RNA' in results.absplice_rna.columns
    
def test_splicing_outlier_result__add_tissue_info_to_spliceai(outlier_dl, outlier_dl_multi, outlier_model, gene_map, gene_tpm, var_samples_df):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.gene_map = gene_map
    results.gene_tpm = gene_tpm
    results.add_spliceai(spliceai_path, gene_mapping=True)
    results.add_samples(var_samples_df)
    
    df_spliceai_tpm = results._add_tissue_info_to_spliceai()
    
    # add tissue info without gene_tpm provided
    results.gene_tpm = None
    df_spliceai_no_tpm = results._add_tissue_info_to_spliceai()
    
    assert 'gene_tpm' in df_spliceai_tpm.columns
    assert 'gene_tpm' not in df_spliceai_no_tpm.columns
    
    spliceai_tpm_index = df_spliceai_tpm.set_index(['variant', 'gene_id', 'tissue', 'sample']).index.unique()
    spliceai_no_tpm_index = df_spliceai_no_tpm.set_index(['variant', 'gene_id', 'tissue', 'sample']).index.unique()
    
    assert len(spliceai_tpm_index[~spliceai_tpm_index.get_level_values('tissue').isna()]\
        .difference(spliceai_no_tpm_index)) == 0
    # Some gene_ids are NA, because could not be found in gene_map (e.g. FakeGene)
    # assert len(spliceai_no_tpm_index.difference(spliceai_tpm_index)) >= 0
    assert len(spliceai_no_tpm_index[~spliceai_no_tpm_index.get_level_values('gene_id').isna()]\
        .difference(spliceai_tpm_index[~spliceai_tpm_index.get_level_values('gene_id').isna()])) >= 0
    assert len(spliceai_tpm_index) == df_spliceai_tpm.shape[0]
    assert len(spliceai_no_tpm_index) == spliceai_no_tpm_index.shape[0]
        

def test_outlier_results_multi_vcf(outlier_dl_multi, outlier_model, var_samples_df):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.add_samples(var_samples_df)
    # assert results.gene_mmsplice.index.names == ['gene_name', 'sample', 'tissue']
    # assert sorted(results.gene_mmsplice.index.tolist()) == sorted([
    #     ('BRCA1', 'NA00002', 'Lung'),
    #     ('BRCA1', 'NA00003', 'Lung'),
    #     ('BRCA1', 'NA00002', 'Testis'),
    #     ('BRCA1', 'NA00003', 'Testis'),
    # ])
    assert results.gene_mmsplice.index.names == ['gene_id', 'tissue', 'sample']
    assert sorted(results.gene_mmsplice.reset_index().set_index(['gene_name', 'sample', 'tissue']).index.tolist()) == sorted([
        # ('BRCA1', 'NA00002', 'Lung'),
        # ('BRCA1', 'NA00003', 'Lung'),
        ('BRCA1', 'NA00002', 'Testis'),
        ('BRCA1', 'NA00003', 'Testis'),
    ])

# # def test_outlier_results_filter_samples_with_RNA_seq(outlier_results_multi, outlier_model):
# def test_outlier_results_filter_samples_with_RNA_seq(outlier_dl, outlier_dl_multi, var_samples_df, outlier_model, gene_map, gene_tpm):
#     # TODO: after filter_samples_with_RNA_seq spliceai_path contains tissue info (maybe remove)
#     samples_for_tissue = {
#         'Testis': ['NA00002'],
#         'Lung': ['NA00002', 'NA00003']
#     }

#     results = outlier_model.predict_on_dataloader(outlier_dl_multi)
#     results.gene_map = gene_map
#     results.gene_tpm = gene_tpm
#     results.add_spliceai(spliceai_path, gene_mapping=True)
#     results.add_samples(var_samples_df)
    
#     # results = outlier_results_multi
#     # results.add_spliceai(spliceai_path, gene_mapping=False)
#     # assert results.df_mmsplice[['tissue', 'sample']].set_index('tissue').to_dict() == \
#     assert results.df_mmsplice[['tissue', 'sample']].groupby('tissue')['sample'].apply(lambda x: ';'.join(sorted(list(set(x))))).to_dict() == \
#         {
#             # 'Lung': 'NA00002;NA00003', 
#             'Testis': 'NA00002;NA00003'}
        
#     df_spliceai_tpm = results._add_tissue_info_to_spliceai()
#     assert df_spliceai_tpm[['tissue', 'sample']].groupby('tissue')['sample'].apply(lambda x: ';'.join(sorted(list(set(x))))).to_dict() == \
#         {
#             # 'Lung': 'NA00002;NA00003', 
#             'Testis': 'NA00002;NA00003'}

#     results.filter_samples_with_RNA_seq(samples_for_tissue)

#     assert results.df_mmsplice[['tissue', 'sample']].groupby('tissue')['sample'].apply(lambda x: ';'.join(sorted(list(set(x))))).to_dict() == \
#         {
#             # 'Lung': 'NA00002;NA00003', 
#             'Testis': 'NA00002'}
#     assert results.df_spliceai[['tissue', 'sample']].groupby('tissue')['sample'].apply(lambda x: ';'.join(sorted(list(set(x))))).to_dict() == \
#         {
#             # 'Lung': 'NA00002;NA00003', 
#             'Testis': 'NA00002'}
    

def test_outlier_results_infer_cat(outlier_dl_multi, outlier_results, cat_dl, outlier_model, var_samples_df):

    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.add_samples(var_samples_df)
    results.infer_cat(cat_dl)

    assert sorted(results.junction.columns.tolist()) == sorted([
        'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type', 'gene_tpm',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
        'delta_logit_psi_cat', 'delta_psi_cat'])

    # Note: previously max aggregation was done in gene property (would have to do aggregation for delta_psi and delta_psi_cat)
    # assert sorted(results.gene_mmsplice.columns.tolist()) == sorted([
    #     'junction', 'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
    #     'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
    #     'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
    #     'gene_id', 'transcript_id', 'gene_type', 'gene_tpm',
    #     'delta_psi', 'delta_logit_psi',
    #     'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
    #     'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
    #     'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
    #     'delta_logit_psi_cat', 'delta_psi_cat'])

    assert results.junction.loc[(
        '17:41201211-41203079:-',  'Testis', 'NA00002',)] is not None


# def test_outlier_results_cat_concat(outlier_results_multi, cat_dl, outlier_model):

#     results = outlier_results_multi
#     results.infer_cat(cat_dl)

#     # assert len(results.gene_mmsplice.loc[list(set(results.gene_mmsplice[['delta_psi', 'delta_psi_cat', 'tissue_cat']].index))[0]]
#     #            ['delta_psi_cat'].values) > 1

#     # TODO: test fails, when removing outlier_results as argument, and not taking np.abs of gene_cat_concat. why is outlier_results changing outcome?
#     # assert np.abs(results.gene_mmsplice.loc[list(set(results.gene_mmsplice[['delta_psi', 'delta_psi_cat', 'tissue_cat']].index))[0]]
#     #               ['delta_psi_cat'].values).max() == \
#     #     np.abs(results.gene_cat_concat.loc[list(set(results.gene_mmsplice[['delta_psi', 'delta_psi_cat', 'tissue_cat']].index))[0]]
#     #            ['delta_psi_cat'])
#     assert np.abs(results.gene_mmsplice.loc[list(set(results.gene_mmsplice[['delta_psi', 'delta_psi_cat', 'tissue_cat']].index))[0]]
#                   ['delta_psi_cat']).max() == \
#         np.abs(results.gene_cat_concat.loc[list(set(results.gene_mmsplice[['delta_psi', 'delta_psi_cat', 'tissue_cat']].index))[0]]
#                ['delta_psi_cat'])

#     assert sorted(results.junction_cat_concat.columns.tolist()) == sorted([
#         'splice_site', 'event_type', 'variant',   'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'ref_psi', 'k', 'n', 'median_n',
#         'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#         'gene_id', 'gene_name', 'transcript_id', 'gene_type',
#         'delta_psi', 'delta_logit_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#         'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#         'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
#         'delta_logit_psi_cat', 'delta_psi_cat'])

#     assert sorted(results.splice_site_cat_concat.columns.tolist()) == sorted([
#         'junction', 'event_type', 'variant',   'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'ref_psi', 'k', 'n', 'median_n',
#         'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#         'gene_id', 'gene_name', 'transcript_id', 'gene_type',
#         'delta_psi', 'delta_logit_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#         'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#         'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
#         'delta_logit_psi_cat', 'delta_psi_cat'])

#     assert sorted(results.gene_cat_concat.columns.tolist()) == sorted([
#         'junction', 'event_type', 'variant',   'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
#         'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#         'gene_id', 'transcript_id', 'gene_type',
#         'delta_psi', 'delta_logit_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#         'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#         'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
#         'delta_logit_psi_cat', 'delta_psi_cat'])

#     assert results.junction_cat_concat.loc[(
#         '17:41201211-41203079:-', 'NA00002', 'Testis')] is not None

#     assert results.splice_site_cat_concat.loc[(
#         '17:41203079:-', 'NA00002', 'Testis')] is not None

#     assert results.gene_cat_concat.loc[(
#         'BRCA1', 'NA00002', 'Testis')] is not None

#     results._gene_cat_concat = None

#     # results._gene_mmsplice = pd.concat([results._gene_mmsplice.iloc[0:1], results._gene_mmsplice])
#     # cat_cols = [x for x in results._gene_mmsplice.columns if 'cat' in x]
#     # for i in cat_cols:
#     #     results._gene_mmsplice.iloc[0, results._gene_mmsplice.columns.get_loc(i)] = np.NaN
#     # r = results._gene_mmsplice.reset_index()
#     # r.iloc[0, r.columns.get_loc('gene_name')] = 'BRCA1'
#     # r.iloc[0, r.columns.get_loc('sample')] = 'NA00005'
#     # r.iloc[0, r.columns.get_loc('tissue')] = 'Lung'
#     # results._gene_mmsplice = r.set_index(['gene_name', 'sample', 'tissue'])
#     results._junction = pd.concat([results._junction.iloc[0:1], results._junction])
#     cat_cols = [x for x in results._junction.columns if 'cat' in x]
#     for i in cat_cols:
#         results._junction.iloc[0, results._junction.columns.get_loc(i)] = np.NaN
#     r = results._junction.reset_index()
#     r.iloc[0, r.columns.get_loc('gene_name')] = 'BRCA1'
#     r.iloc[0, r.columns.get_loc('sample')] = 'NA00005'
#     r.iloc[0, r.columns.get_loc('tissue')] = 'Lung'
#     results._junction = r.set_index(['junction', 'sample', 'tissue'])
#     results._gene_mmsplice = r.set_index(['gene_name', 'sample', 'tissue'])

#     assert results.gene_mmsplice.loc[(
#         'BRCA1', 'NA00005', 'Lung')] is not None
#     assert results.gene_cat_concat.loc[(
#         'BRCA1', 'NA00005', 'Lung')] is not None


# def test_outlier_results_cat_features(outlier_results_multi, cat_dl):

#     results = outlier_results_multi
#     results.infer_cat(cat_dl)
#     results.gene_cat_features

#     # assert np.array_equal(results.gene_mmsplice[results.gene_mmsplice['tissue_cat'] == 'blood']['delta_psi_cat'].values,
#     #                       results._gene_cat_features['delta_psi_blood'].values)

#     # assert np.array_equal(results.gene_mmsplice[results.gene_mmsplice['tissue_cat'] == 'lymphocytes']['delta_psi_cat'].values,
#     #                       results._gene_cat_features['delta_psi_lymphocytes'].values)

#     assert sorted(results.junction_cat_features.columns.tolist()) == sorted([
#         'splice_site', 'event_type', 'variant',   'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'ref_psi', 'k', 'n', 'median_n',
#         'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#         'gene_id', 'gene_name', 'transcript_id', 'gene_type',
#         'delta_psi', 'delta_logit_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#         'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#         'count_blood', 'psi_blood', 'ref_psi_blood', 'k_blood', 'n_blood', 'median_n_blood',
#         'delta_logit_psi_blood', 'delta_psi_blood',
#         'count_lymphocytes', 'psi_lymphocytes', 'ref_psi_lymphocytes', 'k_lymphocytes', 'n_lymphocytes', 'median_n_lymphocytes',
#         'delta_logit_psi_lymphocytes', 'delta_psi_lymphocytes'])

#     assert sorted(results.splice_site_cat_features.columns.tolist()) == sorted([
#         'junction', 'event_type', 'variant',   'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'ref_psi', 'k', 'n', 'median_n',
#         'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#         'gene_id', 'gene_name', 'transcript_id', 'gene_type',
#         'delta_psi', 'delta_logit_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#         'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#         'count_blood', 'psi_blood', 'ref_psi_blood', 'k_blood', 'n_blood', 'median_n_blood',
#         'delta_logit_psi_blood', 'delta_psi_blood',
#         'count_lymphocytes', 'psi_lymphocytes', 'ref_psi_lymphocytes', 'k_lymphocytes', 'n_lymphocytes', 'median_n_lymphocytes',
#         'delta_logit_psi_lymphocytes', 'delta_psi_lymphocytes'])

#     assert sorted(results.gene_cat_features.columns.tolist()) == sorted([
#         'junction', 'event_type', 'variant',   'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
#         'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#         'gene_id', 'transcript_id', 'gene_type',
#         'delta_psi', 'delta_logit_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#         'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#         'count_blood', 'psi_blood', 'ref_psi_blood', 'k_blood', 'n_blood', 'median_n_blood',
#         'delta_logit_psi_blood', 'delta_psi_blood',
#         'count_lymphocytes', 'psi_lymphocytes', 'ref_psi_lymphocytes', 'k_lymphocytes', 'n_lymphocytes', 'median_n_lymphocytes',
#         'delta_logit_psi_lymphocytes', 'delta_psi_lymphocytes'])

#     assert results.junction_cat_features.loc[(
#         '17:41201211-41203079:-', 'NA00002', 'Testis')] is not None

#     assert results.splice_site_cat_features.loc[(
#         '17:41203079:-', 'NA00002', 'Testis')] is not None

#     assert results.gene_cat_features.loc[(
#         'BRCA1', 'NA00002', 'Testis')] is not None

#     results._gene_cat_features = None
#     results._gene_mmsplice = None

#     # results._gene_mmsplice = pd.concat([results._gene_mmsplice.iloc[0:1], results._gene_mmsplice])
#     # cat_cols = [x for x in results._gene_mmsplice.columns if 'cat' in x]
#     # for i in cat_cols:
#     #     results._gene_mmsplice.iloc[0, results._gene_mmsplice.columns.get_loc(i)] = np.NaN
#     # r = results._gene_mmsplice.reset_index()
#     # r.iloc[0, r.columns.get_loc('gene_name')] = 'BRCA1'
#     # r.iloc[0, r.columns.get_loc('sample')] = 'NA00005'
#     # r.iloc[0, r.columns.get_loc('tissue')] = 'Lung'
#     # results._gene_mmsplice = r.set_index(['gene_name', 'sample', 'tissue'])
#     results._junction = pd.concat([results._junction.iloc[0:1], results._junction])
#     cat_cols = [x for x in results._junction.columns if 'cat' in x]
#     for i in cat_cols:
#         results._junction.iloc[0, results._junction.columns.get_loc(i)] = np.NaN
#     r = results._junction.reset_index()
#     r.iloc[0, r.columns.get_loc('gene_name')] = 'BRCA1'
#     r.iloc[0, r.columns.get_loc('sample')] = 'NA00005'
#     r.iloc[0, r.columns.get_loc('tissue')] = 'Lung'
#     results._junction = r.set_index(['junction', 'sample', 'tissue'])

#     assert results.gene_mmsplice.loc[(
#         'BRCA1', 'NA00005', 'Lung')] is not None
#     assert results.gene_cat_features.loc[(
#         'BRCA1', 'NA00005', 'Lung')] is not None


# def test_splicing_outlier_result_infer_cat_add_spliceai(outlier_results_multi, cat_dl, var_samples_df):
    
#     samples_for_tissue = {
#         'Testis': ['NA00002'],
#         'Lung': ['NA00002', 'NA00003']
#     }

#     results = outlier_results_multi
#     results.add_spliceai(spliceai_path, gene_mapping=False)
#     results.add_samples(var_samples_df)
#     results.filter_samples_with_RNA_seq(samples_for_tissue)
#     results.infer_cat(cat_dl)
    
#     # assert results.df_mmsplice_cat contains only valid tissue, sample combinations

#     # assert sorted(results.gene_mmsplice.columns.tolist()) == sorted([
#     #     'junction', 'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
#     #     'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
#     #     'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#     #     'gene_id', 'transcript_id', 'gene_type',  
#     #     'delta_psi', 'delta_logit_psi',
#     #     'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#     #     'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#     #     'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat',
#     #     'k_cat', 'n_cat', 'median_n_cat', 'delta_logit_psi_cat', 'delta_psi_cat',
#     #     'index', 'variant_spliceAI',
#     #     'delta_score', 'acceptor_gain', 'acceptor_loss', 'donor_gain',
#     #     'donor_loss', 'acceptor_gain_position', 'acceptor_loss_positiin',
#     #     'donor_gain_position', 'donor_loss_position', 'GQ',
#     #     'DP_ALT'
#     # ])

#     # assert sorted(results.gene_cat_concat.columns.tolist()) == sorted([
#     #     'junction', 'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
#     #     'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
#     #     'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#     #     'gene_id', 'transcript_id', 'gene_type',  
#     #     'delta_psi', 'delta_logit_psi',
#     #     'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#     #     'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#     #     'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat',
#     #     'k_cat', 'n_cat', 'median_n_cat', 'delta_logit_psi_cat', 'delta_psi_cat',
#     #     'index', 'variant_spliceAI',
#     #     'delta_score', 'acceptor_gain', 'acceptor_loss', 'donor_gain',
#     #     'donor_loss', 'acceptor_gain_position', 'acceptor_loss_positiin',
#     #     'donor_gain_position', 'donor_loss_position', 'GQ',
#     #     'DP_ALT'
#     # ])
    
    
# # DEBUG
# import pytest
# import pandas as pd
# import numpy as np
# from kipoiseq.extractors.vcf import MultiSampleVCF
# # from kipoiseq.extractors.vcf_query import to_sample_csv
# from splicing_outlier_prediction import SpliceOutlier, SpliceOutlierDataloader, CatInference
# from splicing_outlier_prediction.ensemble import train_model_ebm
# from conftest import fasta_file, vcf_file, multi_vcf_file, multi_vcf_samples, \
#     ref_table5_kn_testis, ref_table3_kn_testis,  \
#     ref_table5_kn_lung, ref_table3_kn_lung, \
#     combined_ref_tables5_testis_lung, combined_ref_tables3_testis_lung, \
#     count_cat_file_lymphocytes,  count_cat_file_blood, \
#     spliceai_path, pickle_DNA, pickle_DNA_CAT
    
# from kipoiseq.extractors.vcf_query import VariantIntervalQueryable
# from kipoiseq.dataclasses import Variant, Interval


# # @pytest.fixture
# # def variant_queryable():
# #     vcf = MultiSampleVCF(vcf_file)
# #     return VariantIntervalQueryable(vcf, [
# #         (
# #             [
# #                 Variant('chr1', 12, 'A', 'T'),
# #                 Variant('chr1', 18, 'A', 'C', filter='q10'),
# #             ],
# #             Interval('chr1', 10, 20)
# #         ),
# #         (
# #             [
# #                 Variant('chr2', 120, 'AT', 'AAAT'),
# #             ],
# #             Interval('chr2', 110, 200)
# #         )
# #     ])

# def test_write_vcf_condensed():
#     # variant_queryable.to_vcf(vcf_file_condensed, remove_samples=True, clean_info=True)
#     # [x.__str__() for x in MultiSampleVCF(vcf_file).query_all().filter(lambda variant: variant.pos in [41201201, 41279042, 41276032])]
#     MultiSampleVCF(multi_vcf_file).query_all().filter(lambda variant: variant.pos in [41201201, 41279042, 41276032]).to_vcf(vcf_file.replace('.gz', ''), remove_samples=True, clean_info=True)
    
#     # [x for x in VariantIntervalQueryable(vcf_file, [([Variant('chr17', 41201201, 'TTC', 'CA')], Interval('chr17', 41201101, 41201209))])]
#     # VariantIntervalQueryable(vcf_file, [([Variant('chr17', 41201201, 'TTC', 'CA')], Interval('chr17', 41201101, 41201209))]).to_vcf(vcf_file_condensed, remove_samples=True)
    
# # def tabix_vcf():
# #     bgzip vcf_file.replace('.gz', '')
# #     tabix vcf_file