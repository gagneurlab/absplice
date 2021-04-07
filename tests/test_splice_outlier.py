import pytest
import pandas as pd
from splicing_outlier_prediction import SpliceOutlier, SpliceOutlierDataloader
from splicing_outlier_prediction.spliceAI import SpliceAI
from conftest import ref_table5_kn_file, ref_table3_kn_file, fasta_file, \
    multi_vcf_file, spliceai_db_path


@pytest.fixture
def outlier_model():
    return SpliceOutlier()


def test_splicing_outlier_on_batch(outlier_model, outlier_dl):
    batch = next(outlier_dl.batch_iter())
    df = outlier_model.predict_on_batch(batch, outlier_dl)
    assert df.columns.tolist() == [
        'variant', 'junction', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n',
        'gene_id', 'gene_name', 'weak',
        'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
        'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
        'alt_exon', 'alt_donor', 'alt_donorIntron']
    assert df.shape[0] == 32


def test_splicing_outlier_predict_on_dataloader(outlier_model, outlier_dl):
    results = outlier_model.predict_on_dataloader(outlier_dl)
    results.df.columns.tolist() == [
        'variant', 'junction', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n',
        'gene_id', 'gene_name', 'weak',
        'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
        'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
        'alt_exon', 'alt_donor', 'alt_donorIntron']


def test_splicing_outlier_predict_save(outlier_model, outlier_dl, tmp_path):
    output_csv = tmp_path / 'pred.csv'
    outlier_model.predict_save(outlier_dl, output_csv)
    df = pd.read_csv(output_csv)
    assert df.columns.tolist() == [
        'variant', 'junction', 'event_type',
        'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n',
        'gene_id', 'gene_name', 'weak',
        'transcript_id', 'gene_type',
        'delta_logit_psi', 'delta_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
        'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
        'alt_exon', 'alt_donor', 'alt_donorIntron']


@pytest.fixture
def outlier_results(outlier_model, outlier_dl):
    return outlier_model.predict_on_dataloader(outlier_dl)


def test_splicing_outlier_result_psi(outlier_results):
    results = outlier_results.psi5
    assert all(results.psi5.df['event_type'] == 'psi5')

    results = outlier_results.psi3
    assert all(results.psi3.df['event_type'] == 'psi3')


def test_splicing_outlier_result_splice_site(outlier_results):
    assert sorted(outlier_results.splice_site.index) == sorted(
        set(outlier_results.df['splice_site']))


def test_splicing_outlier_result_gene(outlier_results):
    assert sorted(outlier_results.gene.index) == sorted(
        set(outlier_results.df['gene_name']))


def test_outlier_results_multi_vcf(outlier_model):
    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_file], ref_tables3=[ref_table3_kn_file],
        samples=True)

    results = outlier_model.predict_on_dataloader(dl)
    assert results.gene.index.names == ['gene_name', 'sample']
    assert results.gene.index.tolist() == [
        ('BRCA1', 'NA00002'),
        ('BRCA1', 'NA00003')
    ]


def test_outlier_results_infer_cat_add_spliceAI(outlier_results, cat_dl, outlier_model):
    with pytest.raises(ValueError):
        outlier_results.infer_cat(cat_dl)

    dl = SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_file], ref_tables3=[ref_table3_kn_file],
        samples=True)

    results = outlier_model.predict_on_dataloader(dl)

    spliceai = SpliceAI(db_path=spliceai_db_path)
    df_spliceAI = spliceai.predict_df(
        ['17:34149615:A>T', '17:34149617:A>T']).reset_index()
    df_spliceAI.loc[0, 'samples'] = 'NA00001;NA00002;NA00003'
    results.add_spliceAI(df_spliceAI)
    results.df.loc[0, 'samples'] = 'NA00002;NA00003;NA00004'
    results.infer_cat(cat_dl)

    assert sorted(results.junction.columns.tolist()) == sorted([
        'event_type', 'variant', 'genotype', 'Chromosome', 'Start', 'End', 'Strand',
        'events', 'splice_site', 'ref_psi', 'k', 'n', 'gene_id', 'gene_name',
        'weak', 'transcript_id', 'gene_type', 'delta_psi', 'delta_logit_psi',
        'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor',
        'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor', 'alt_exon',
        'alt_donor', 'alt_donorIntron', 'count_cat', 'psi_cat', 'ref_psi_cat',
        'k_cat', 'n_cat', 'delta_logit_psi_cat', 'delta_psi_cat'])

    assert sorted(results.gene.columns.tolist()) == sorted([
        'junction', 'event_type', 'variant', 'genotype', 'Chromosome', 'Start',
        'End', 'Strand', 'events', 'splice_site', 'ref_psi', 'k', 'n',
        'gene_id', 'weak', 'transcript_id', 'delta_psi', 'gene_type',
        'delta_logit_psi', 'ref_acceptorIntron', 'ref_acceptor', 'ref_exon',
        'ref_donor', 'ref_donorIntron', 'alt_acceptorIntron', 'alt_acceptor',
        'alt_exon', 'alt_donor', 'alt_donorIntron', 'count_cat', 'psi_cat',
        'ref_psi_cat', 'k_cat', 'n_cat', 'delta_logit_psi_cat', 'delta_psi_cat',
        'index', 'variant_spliceAI', 'delta_score', 'acceptor_gain',
        'acceptor_loss', 'donor_gain', 'donor_loss', 'acceptor_gain_position',
        'acceptor_loss_positiin', 'donor_gain_position', 'donor_loss_position'])

    assert results.junction.loc[(
        '17:41201211-41203079:-', 'NA00004')] is not None
