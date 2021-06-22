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


def test_splicing_outlier_predict_on_dataloader_correct_tissue(outlier_model, outlier_dl):
    results = outlier_model.predict_on_dataloader(outlier_dl)
    assert False not in results.df.apply(lambda x: x['event_type'] in x['tissue'], axis=1)


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
