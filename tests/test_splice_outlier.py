import pytest
import pandas as pd
from splicing_outlier_prediction import SpliceOutlier, SpliceOutlierDataloader, CatInference
from conftest import ref_table5_kn_file, ref_table3_kn_file, fasta_file, \
    multi_vcf_file, spliceai_db_path,count_cat_file, ref_table5_kn_file2, \
    ref_table3_kn_file2, combined_ref_tables5_file, combined_ref_tables3_file, \
    ref_table5_kn_file_fake1, ref_table5_kn_file_fake2, \
    ref_table3_kn_file_fake1, ref_table3_kn_file_fake2, \
    combined_ref_tables5_file_fake, combined_ref_tables3_file_fake


# @pytest.fixture
# def outlier_model():
#     return SpliceOutlier()


def test_splicing_outlier_on_batch(outlier_model, outlier_dl):
    batch = next(outlier_dl.batch_iter())
    df = outlier_model.predict_on_batch(batch, outlier_dl)
    assert sorted(df.columns.tolist()) == sorted([
        'variant', 'junction', 'tissue', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n',
        'weak_site_acceptor', 'median_n', 'novel_junction', 'weak_site_donor',
        'gene_id', 'gene_name', 'weak',
        'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
        'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
        'alt_exon', 'alt_donor', 'alt_donorIntron'])
    assert df.shape[0] == 32


def test_splicing_outlier_predict_on_dataloader(outlier_model, outlier_dl):
    results = outlier_model.predict_on_dataloader(outlier_dl)
    assert sorted(results.df.columns.tolist()) == sorted([
        'variant', 'junction', 'tissue', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n',
        'weak_site_acceptor', 'median_n', 'novel_junction', 'weak_site_donor',
        'gene_id', 'gene_name', 'weak',
        'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
        'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
        'alt_exon', 'alt_donor', 'alt_donorIntron'])


def test_splicing_outlier_predict_save(outlier_model, outlier_dl, tmp_path):
    output_csv = tmp_path / 'pred.csv'
    outlier_model.predict_save(outlier_dl, output_csv)
    df = pd.read_csv(output_csv)
    assert sorted(df.columns.tolist()) == sorted([
        'variant', 'junction', 'tissue', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 
        'weak_site_acceptor', 'median_n', 'novel_junction', 'weak_site_donor',
        'gene_id', 'gene_name', 'weak',
        'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
        'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
        'alt_exon', 'alt_donor', 'alt_donorIntron'])


# @pytest.fixture
# def outlier_results(outlier_model, outlier_dl):
#     return outlier_model.predict_on_dataloader(outlier_dl)


def test_tissues(outlier_results):
    assert sorted(outlier_results.df['tissue'].unique()) == sorted(['lymphocytes', 'lung'])


def test_splicing_outlier_result_psi(outlier_results):
    results = outlier_results.psi5
    assert all(results.psi5.df['event_type'] == 'psi5')

    results = outlier_results.psi3
    assert all(results.psi3.df['event_type'] == 'psi3')


def test_splicing_outlier_result_splice_site(outlier_results):
    assert sorted(outlier_results.splice_site.index.get_level_values('splice_site')) \
        == sorted(set(outlier_results.df['splice_site']))


def test_splicing_outlier_result_gene(outlier_results):
    assert sorted(outlier_results.gene.index.get_level_values('gene_name')) \
        == sorted(set(outlier_results.df['gene_name']))


def test_outlier_results_multi_vcf(outlier_model):
    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_file], 
        ref_tables3=[ref_table3_kn_file],
        samples=True,
        regex_pattern='test_(.*)_ref'
    )

    results = outlier_model.predict_on_dataloader(dl)
    assert results.gene.index.names == ['gene_name', 'sample', 'tissue']
    assert results.gene.index.tolist() == [
        ('BRCA1', 'NA00002', 'lymphocytes'),
        ('BRCA1', 'NA00003', 'lymphocytes')
    ]


def test_outlier_results_infer_cat(outlier_results, cat_dl, outlier_model):
# def test_outlier_results_infer_cat(outlier_model):

    # with pytest.raises(ValueError):
    #     outlier_results.infer_cat(cat_dl)

    # cat_dl1 = CatInference(
    #    # ref_tables5=[ref_table5_kn_file, ref_table5_kn_file2], 
    #    # ref_tables3=[ref_table3_kn_file, ref_table3_kn_file2],
    #     ref_tables5=[ref_table5_kn_file_fake1, ref_table5_kn_file_fake2], 
    #     ref_tables3=[ref_table3_kn_file_fake1, ref_table3_kn_file_fake2],
    #     regex_pattern='test_(.*)_ref',
    #     count_cat=count_cat_file)

    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_file, ref_table5_kn_file2], 
        ref_tables3=[ref_table3_kn_file, ref_table3_kn_file2],
        combined_ref_tables5=combined_ref_tables5_file, 
        combined_ref_tables3=combined_ref_tables3_file,
        # ref_tables5=[ref_table5_kn_file_fake1, ref_table5_kn_file_fake2], 
        # ref_tables3=[ref_table3_kn_file_fake1, ref_table3_kn_file_fake2],
        # combined_ref_tables5=combined_ref_tables5_file_fake, 
        # combined_ref_tables3=combined_ref_tables3_file_fake,
        regex_pattern='test_(.*)_ref',
        samples=True)

    results = outlier_model.predict_on_dataloader(dl)
    results.infer_cat(cat_dl)

    assert sorted(results.junction.columns.tolist()) == sorted([
        'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'gene_id', 'gene_name',
        'weak', 'transcript_id', 'gene_type', 'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
        'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor', 'alt_exon',
        'alt_donor', 'alt_donorIntron', 'count_cat', 'psi_cat', 'ref_psi_cat',
        'k_cat', 'n_cat', 'median_n_cat', 'delta_logit_psi_cat', 'delta_psi_cat'])

    assert sorted(results.gene.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start',
        'End', 'Strand', 'events', 'splice_site', 'ref_psi', 'k', 'n',
        'gene_id', 'weak', 'transcript_id', 'delta_psi', 'gene_type',
        'delta_logit_psi', 'ref_acceptorIntron', 'ref_acceptor', 'ref_exon',
        'ref_donor', 'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
        'alt_exon', 'alt_donor', 'alt_donorIntron', 'count_cat', 'psi_cat',
        'ref_psi_cat', 'k_cat', 'n_cat', 'median_n_cat', 'delta_logit_psi_cat', 'delta_psi_cat',
        'index', 'variant_spliceAI', 'delta_score', 'acceptor_gain',
        'acceptor_loss', 'donor_gain', 'donor_loss', 'acceptor_gain_position',
        'acceptor_loss_positiin', 'donor_gain_position', 'donor_loss_position'])

    assert results.junction.loc[(
        '17:41201211-41203079:-', 'NA00004')] is not None
