import pytest
import pandas as pd
import numpy as np
from splicing_outlier_prediction import SpliceOutlier, SpliceOutlierDataloader, CatInference
from splicing_outlier_prediction.ensemble import train_model_ebm
from conftest import fasta_file, multi_vcf_file, \
    ref_table5_kn_testis, ref_table3_kn_testis,  \
    ref_table5_kn_lung, ref_table3_kn_lung, \
    combined_ref_tables5_testis_lung, combined_ref_tables3_testis_lung, \
    count_cat_file_lymphocytes,  count_cat_file_blood, \
    spliceAI, pickle_DNA, pickle_DNA_CAT


def test_tissues(outlier_results):
    assert sorted(outlier_results.df['tissue'].unique()) == sorted(
        ['testis', 'lung'])


def test_splicing_outlier_result_psi(outlier_results):
    results = outlier_results.psi5
    assert all(results.psi5.df['event_type'] == 'psi5')

    results = outlier_results.psi3
    assert all(results.psi3.df['event_type'] == 'psi3')


def test_splicing_outlier_result_splice_site(outlier_results):
    assert sorted(outlier_results.splice_site.index) \
        == sorted(set(outlier_results.df.set_index(['splice_site', 'tissue']).index))


def test_splicing_outlier_result_gene(outlier_results):
    assert sorted(set(outlier_results.gene.index)) \
        == sorted(set(outlier_results.df.set_index(['gene_name', 'tissue']).index))


def test_outlier_results_multi_vcf(outlier_model):
    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_testis],
        ref_tables3=[ref_table3_kn_testis],
        samples=True,
        regex_pattern='test_(.*)_ref'
    )

    results = outlier_model.predict_on_dataloader(dl)
    assert results.gene.index.names == ['gene_name', 'sample', 'tissue']
    assert results.gene.index.tolist() == [
        ('BRCA1', 'NA00002', 'testis'),
        ('BRCA1', 'NA00003', 'testis')
    ]


def test_outlier_results_filter_samples_with_RNA_seq(outlier_model):
    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung],
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        combined_ref_tables5=combined_ref_tables5_testis_lung,
        combined_ref_tables3=combined_ref_tables3_testis_lung,
        regex_pattern='test_(.*)_ref',
        samples=True)

    samples_for_tissue = {
        'testis': ['NA00002'],
        'lung': ['NA00002', 'NA00003']
    }

    results = outlier_model.predict_on_dataloader(dl)
    results.add_spliceAI(spliceAI)
    assert results.df[['tissue', 'samples']].set_index('tissue').to_dict() == \
        {'samples': {'lung': 'NA00002;NA00003', 'testis': 'NA00002;NA00003'}}
    assert results.df_spliceAI[results.df_spliceAI['variant'] == '17:41154337:ACT>A']['samples'].values == \
        'NA00002;NA00003;NA00001'
    assert results.df_spliceAI.shape[0] == 34

    results.filter_samples_with_RNA_seq(samples_for_tissue)

    assert results.df[['tissue', 'samples']].set_index('tissue').to_dict() == \
        {'samples': {'lung': 'NA00002;NA00003', 'testis': 'NA00002'}}
    assert results.df_spliceAI[results.df_spliceAI['variant'] == '17:41154337:ACT>A'][['tissue', 'samples']].set_index('tissue').to_dict() == \
        {'samples': {'lung': 'NA00002;NA00003', 'testis': 'NA00002'}}
    # not twice the size, because some variants only exist for missing samples in target tissues
    assert results.df_spliceAI.shape[0] == 57


def test_outlier_results_infer_cat(outlier_results, cat_dl, outlier_model):
    # with pytest.raises(ValueError):
    #     outlier_results.infer_cat(cat_dl)

    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung],
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        combined_ref_tables5=combined_ref_tables5_testis_lung,
        combined_ref_tables3=combined_ref_tables3_testis_lung,
        regex_pattern='test_(.*)_ref',
        samples=True)

    results = outlier_model.predict_on_dataloader(dl)
    results.infer_cat(cat_dl)

    assert sorted(results.junction.columns.tolist()) == sorted([
        'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
        'delta_logit_psi_cat', 'delta_psi_cat'])

    assert sorted(results.gene.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'transcript_id', 'gene_type',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
        'delta_logit_psi_cat', 'delta_psi_cat'])
    assert results.junction.loc[(
        '17:41201211-41203079:-', 'NA00002', 'testis')] is not None


def test_outlier_results_cat_concat(cat_dl, outlier_model):
    # with pytest.raises(ValueError):
    #     outlier_results.infer_cat(cat_dl)

    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung],
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        combined_ref_tables5=combined_ref_tables5_testis_lung,
        combined_ref_tables3=combined_ref_tables3_testis_lung,
        save_combined_ref_tables=False,
        regex_pattern='test_(.*)_ref',
        samples=True)

    results = outlier_model.predict_on_dataloader(dl)
    results.infer_cat(cat_dl)

    assert len(results.gene.loc[list(set(results.gene[['delta_psi', 'delta_psi_cat', 'tissue_cat']].index))[0]]
               ['delta_psi_cat'].values) > 1

    # TODO: test fails, when removing outlier_results as argument, and not taking np.abs of gene_cat_concat. why is outlier_results changing outcome?
    assert np.abs(results.gene.loc[list(set(results.gene[['delta_psi', 'delta_psi_cat', 'tissue_cat']].index))[0]]
                  ['delta_psi_cat'].values).max() == \
        np.abs(results.gene_cat_concat.loc[list(set(results.gene[['delta_psi', 'delta_psi_cat', 'tissue_cat']].index))[0]]
               ['delta_psi_cat'])

    assert sorted(results.junction_cat_concat.columns.tolist()) == sorted([
        'splice_site', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
        'delta_logit_psi_cat', 'delta_psi_cat'])

    assert sorted(results.splice_site_cat_concat.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
        'delta_logit_psi_cat', 'delta_psi_cat'])

    assert sorted(results.gene_cat_concat.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'transcript_id', 'gene_type',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat',
        'delta_logit_psi_cat', 'delta_psi_cat'])

    assert results.junction_cat_concat.loc[(
        '17:41201211-41203079:-', 'NA00002', 'testis')] is not None

    assert results.splice_site_cat_concat.loc[(
        '17:41203079:-', 'NA00002', 'testis')] is not None

    assert results.gene_cat_concat.loc[(
        'BRCA1', 'NA00002', 'testis')] is not None

    results._gene_cat_concat = None

    results._gene = pd.concat([results._gene.iloc[0:1], results._gene])
    cat_cols = [x for x in results._gene.columns if 'cat' in x]
    for i in cat_cols:
        results._gene.iloc[0, results._gene.columns.get_loc(i)] = np.NaN
    r = results._gene.reset_index()
    r.iloc[0, r.columns.get_loc('gene_name')] = 'BRCA1'
    r.iloc[0, r.columns.get_loc('sample')] = 'NA00005'
    r.iloc[0, r.columns.get_loc('tissue')] = 'lung'
    results._gene = r.set_index(['gene_name', 'sample', 'tissue'])

    assert results.gene.loc[(
        'BRCA1', 'NA00005', 'lung')] is not None
    assert results.gene_cat_concat.loc[(
        'BRCA1', 'NA00005', 'lung')] is not None


def test_outlier_results_cat_features(outlier_results, cat_dl, outlier_model):
    # with pytest.raises(ValueError):
    #     outlier_results.infer_cat(cat_dl)

    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung],
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        combined_ref_tables5=combined_ref_tables5_testis_lung,
        combined_ref_tables3=combined_ref_tables3_testis_lung,
        regex_pattern='test_(.*)_ref',
        samples=True)

    results = outlier_model.predict_on_dataloader(dl)
    results.infer_cat(cat_dl)
    results.gene_cat_features

    assert np.array_equal(results.gene[results.gene['tissue_cat'] == 'blood']['delta_psi_cat'].values,
                          results._gene_cat_features['delta_psi_blood'].values)

    assert np.array_equal(results.gene[results.gene['tissue_cat'] == 'lymphocytes']['delta_psi_cat'].values,
                          results._gene_cat_features['delta_psi_lymphocytes'].values)

    assert sorted(results.junction_cat_features.columns.tolist()) == sorted([
        'splice_site', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'count_blood', 'psi_blood', 'ref_psi_blood', 'k_blood', 'n_blood', 'median_n_blood',
        'delta_logit_psi_blood', 'delta_psi_blood',
        'count_lymphocytes', 'psi_lymphocytes', 'ref_psi_lymphocytes', 'k_lymphocytes', 'n_lymphocytes', 'median_n_lymphocytes',
        'delta_logit_psi_lymphocytes', 'delta_psi_lymphocytes'])

    assert sorted(results.splice_site_cat_features.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'count_blood', 'psi_blood', 'ref_psi_blood', 'k_blood', 'n_blood', 'median_n_blood',
        'delta_logit_psi_blood', 'delta_psi_blood',
        'count_lymphocytes', 'psi_lymphocytes', 'ref_psi_lymphocytes', 'k_lymphocytes', 'n_lymphocytes', 'median_n_lymphocytes',
        'delta_logit_psi_lymphocytes', 'delta_psi_lymphocytes'])

    assert sorted(results.gene_cat_features.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'transcript_id', 'gene_type',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'count_blood', 'psi_blood', 'ref_psi_blood', 'k_blood', 'n_blood', 'median_n_blood',
        'delta_logit_psi_blood', 'delta_psi_blood',
        'count_lymphocytes', 'psi_lymphocytes', 'ref_psi_lymphocytes', 'k_lymphocytes', 'n_lymphocytes', 'median_n_lymphocytes',
        'delta_logit_psi_lymphocytes', 'delta_psi_lymphocytes'])

    assert results.junction_cat_features.loc[(
        '17:41201211-41203079:-', 'NA00002', 'testis')] is not None

    assert results.splice_site_cat_features.loc[(
        '17:41203079:-', 'NA00002', 'testis')] is not None

    assert results.gene_cat_features.loc[(
        'BRCA1', 'NA00002', 'testis')] is not None

    results._gene_cat_features = None

    results._gene = pd.concat([results._gene.iloc[0:1], results._gene])
    cat_cols = [x for x in results._gene.columns if 'cat' in x]
    for i in cat_cols:
        results._gene.iloc[0, results._gene.columns.get_loc(i)] = np.NaN
    r = results._gene.reset_index()
    r.iloc[0, r.columns.get_loc('gene_name')] = 'BRCA1'
    r.iloc[0, r.columns.get_loc('sample')] = 'NA00005'
    r.iloc[0, r.columns.get_loc('tissue')] = 'lung'
    results._gene = r.set_index(['gene_name', 'sample', 'tissue'])

    assert results.gene.loc[(
        'BRCA1', 'NA00005', 'lung')] is not None
    assert results.gene_cat_features.loc[(
        'BRCA1', 'NA00005', 'lung')] is not None


def test_splicing_outlier_result_add_spliceAI(outlier_results, outlier_model):
    # with pytest.raises(ValueError):
    #     outlier_results.infer_cat(cat_dl)

    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung],
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        combined_ref_tables5=combined_ref_tables5_testis_lung,
        combined_ref_tables3=combined_ref_tables3_testis_lung,
        regex_pattern='test_(.*)_ref',
        samples=True)

    results = outlier_model.predict_on_dataloader(dl)
    results.add_spliceAI(spliceAI)

    assert sorted(results.gene.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'transcript_id', 'gene_type', 'genotype',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'index', 'variant_spliceAI',
        'delta_score', 'acceptor_gain', 'acceptor_loss', 'donor_gain',
        'donor_loss', 'acceptor_gain_position', 'acceptor_loss_positiin',
        'donor_gain_position', 'donor_loss_position', 'GQ',
        'DP_ALT'
    ])

    assert results.gene['delta_score'] is not None


def test_splicing_outlier_result_infer_cat_add_spliceAI(outlier_results, cat_dl, outlier_model):
    # with pytest.raises(ValueError):
    #     outlier_results.infer_cat(cat_dl)

    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung],
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        combined_ref_tables5=combined_ref_tables5_testis_lung,
        combined_ref_tables3=combined_ref_tables3_testis_lung,
        regex_pattern='test_(.*)_ref',
        samples=True)

    samples_for_tissue = {
        'testis': ['NA00002'],
        'lung': ['NA00002', 'NA00003']
    }

    results = outlier_model.predict_on_dataloader(dl)
    results.add_spliceAI(spliceAI)
    results.filter_samples_with_RNA_seq(samples_for_tissue)
    results.infer_cat(cat_dl)

    assert sorted(results.gene.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'transcript_id', 'gene_type', 'genotype',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat',
        'k_cat', 'n_cat', 'median_n_cat', 'delta_logit_psi_cat', 'delta_psi_cat',
        'index', 'variant_spliceAI',
        'delta_score', 'acceptor_gain', 'acceptor_loss', 'donor_gain',
        'donor_loss', 'acceptor_gain_position', 'acceptor_loss_positiin',
        'donor_gain_position', 'donor_loss_position', 'GQ',
        'DP_ALT'
    ])

    assert sorted(results.gene_cat_concat.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
        'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
        'gene_id', 'transcript_id', 'gene_type', 'genotype',
        'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
        'tissue_cat', 'count_cat', 'psi_cat', 'ref_psi_cat',
        'k_cat', 'n_cat', 'median_n_cat', 'delta_logit_psi_cat', 'delta_psi_cat',
        'index', 'variant_spliceAI',
        'delta_score', 'acceptor_gain', 'acceptor_loss', 'donor_gain',
        'donor_loss', 'acceptor_gain_position', 'acceptor_loss_positiin',
        'donor_gain_position', 'donor_loss_position', 'GQ',
        'DP_ALT'
    ])
