import pytest
import pandas as pd
from splicing_outlier_prediction import SpliceOutlier, SpliceOutlierDataloader, CatInference
from conftest import fasta_file, multi_vcf_file, \
    ref_table5_kn_testis, ref_table3_kn_testis,  \
    ref_table5_kn_lung, ref_table3_kn_lung, \
    combined_ref_tables5_testis_lung, combined_ref_tables3_testis_lung, \
    count_cat_file_lymphocytes,  count_cat_file_blood

# @pytest.fixture
# def outlier_model():
#     return SpliceOutlier()


def test_splicing_outlier_on_batch(outlier_model, outlier_dl):
    batch = next(outlier_dl.batch_iter())
    df = outlier_model.predict_on_batch(batch, outlier_dl)
    assert sorted(df.columns.tolist()) == sorted([
        'variant', 'junction', 'tissue', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n', 
        'novel_junction', 'weak_site_donor', 'weak_site_acceptor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron', 
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron'])
    assert df.shape[0] == 49


def test_splicing_outlier_predict_on_dataloader(outlier_model, outlier_dl):
    results = outlier_model.predict_on_dataloader(outlier_dl)
    assert sorted(results.df.columns.tolist()) == sorted([
        'variant', 'junction', 'tissue', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n', 
        'novel_junction', 'weak_site_donor', 'weak_site_acceptor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron', 
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron'])


def test_splicing_outlier_predict_save(outlier_model, outlier_dl, tmp_path):
    output_csv = tmp_path / 'pred.csv'
    outlier_model.predict_save(outlier_dl, output_csv)
    df = pd.read_csv(output_csv)
    assert sorted(df.columns.tolist()) == sorted([
        'variant', 'junction', 'tissue', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n', 
        'novel_junction', 'weak_site_donor', 'weak_site_acceptor',
        'gene_id', 'gene_name', 'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron', 
        'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron'])


def test_tissues(outlier_results):
    assert sorted(outlier_results.df['tissue'].unique()) == sorted(['testis', 'lung'])


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