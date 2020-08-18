import pytest
import pandas as pd
from splicing_outlier_prediction import CatSpliceOutlier, CatDataloader
from conftest import ref_table5_kn_file, ref_table3_kn_file, fasta_file, multi_vcf_file, count_cat_file


@pytest.fixture
def cat_outlier_model():
    return CatSpliceOutlier()


def test_cat_splicing_outlier_on_batch(cat_outlier_model, cat_dl):
    batch = next(cat_dl.batch_iter())
    df = cat_outlier_model.predict_on_batch(batch, cat_dl)
    assert df.columns.tolist() == [
        'target_psi_ref', 'cat_psi_ref', 'NA00001_cat_psi',
        'NA00001_cat_delta_logit_psi', 'NA00002_cat_psi',
        'NA00002_cat_delta_logit_psi', 'NA00003_cat_psi',
        'NA00003_cat_delta_logit_psi', 'NA00001_target_delta_psi',
        'NA00002_target_delta_psi', 'NA00003_target_delta_psi']
    assert df.shape[0] == 32


# def test_cat_splicing_outlier_predict_on_dataloader(outlier_model, outlier_dl):
#     results = outlier_model.predict_on_dataloader(outlier_dl)
#     results.df.columns.tolist() == [
#         'variant', 'junction', 'event_type',
#         'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'splice_site', 'ref_psi', 'k', 'n',
#         'gene_id', 'gene_name', 'weak',
#         'transcript_id', 'gene_type',
#         'delta_logit_psi', 'delta_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
#         'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
#         'alt_exon', 'alt_donor', 'alt_donorIntron']


# def test_cat_splicing_outlier_predict_save(outlier_model, outlier_dl, tmp_path):
#     output_csv = tmp_path / 'pred.csv'
#     outlier_model.predict_save(outlier_dl, output_csv)
#     df = pd.read_csv(output_csv)
#     assert df.columns.tolist() == [
#         'variant', 'junction', 'event_type',
#         'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'splice_site', 'ref_psi', 'k', 'n',
#         'gene_id', 'gene_name', 'weak',
#         'transcript_id', 'gene_type',
#         'delta_logit_psi', 'delta_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
#         'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
#         'alt_exon', 'alt_donor', 'alt_donorIntron']


# @pytest.fixture
# def outlier_results(outlier_model, outlier_dl):
#     return outlier_model.predict_on_dataloader(outlier_dl)


# def test_splicing_outlier_result_splice_site(outlier_results):
#     assert sorted(outlier_results.splice_site.index) == sorted(
#         set(outlier_results.df['splice_site']))


# def test_splicing_outlier_result_gene(outlier_results):
#     assert sorted(outlier_results.gene.index) == sorted(
#         set(outlier_results.df['gene_id']))


# def test_outlier_results_multi_vcf(outlier_model):
#     dl = SpliceOutlierDataloader(
#         fasta_file, multi_vcf_file,
#         ref_table5=ref_table5_kn_file, ref_table3=ref_table3_kn_file,
#         samples=True)

#     results = outlier_model.predict_on_dataloader(dl)
#     assert results.gene.index.names == ['gene_id', 'sample']
#     assert results.gene.index.tolist() == [
#         ('ENSG00000012048.22_5', 'NA00002'),
#         ('ENSG00000012048.22_5', 'NA00003')
#     ]


# def test_outlier_results_multi_vcf_cat(outlier_model):
#     dl = SpliceOutlierDataloader(
#         fasta_file, multi_vcf_file,
#         ref_table5=ref_table5_kn_file, ref_table3=ref_table3_kn_file, count_cat=count_cat_file,
#         samples=True)

#     results = outlier_model.predict_on_dataloader(dl)
#     assert results.gene.index.names == ['gene_id', 'sample']
#     assert results.gene.index.tolist() == [
#         ('ENSG00000012048.22_5', 'NA00002'),
#         ('ENSG00000012048.22_5', 'NA00003')
#     ]
