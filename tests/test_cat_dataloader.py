import pytest
import numpy as np
from splicing_outlier_prediction import CatInference
from conftest import ref_table5_kn_file, ref_table3_kn_file, fasta_file, \
    vcf_file, multi_vcf_file, count_cat_file


@pytest.fixture
def cat_dl5():
    return CatInference(
        ref_table5=ref_table5_kn_file,
        count_cat=count_cat_file)


@pytest.fixture
def cat_dl3():
    return CatInference(
        ref_table3=ref_table3_kn_file,
        count_cat=count_cat_file)


def test_cat_dataloader_init(cat_dl):
    assert cat_dl.ref_table5.method == 'kn'
    assert cat_dl.ref_table5.df.shape[0] == 26

    assert cat_dl.ref_table3.method == 'kn'
    assert cat_dl.ref_table3.df.shape[0] == 28


def test_cat_dataloader_common(cat_dl):
    assert len(cat_dl.common_junctions5) == 25
    assert len(cat_dl.common_junctions3) == 27


def test_cat_dataloader_common5(cat_dl5):
    assert len(cat_dl5.common_junctions5) == 25


def test_cat_dataloader_common3(cat_dl3):
    assert len(cat_dl3.common_junctions3) == 27


def test_cat_dataloader_infer(cat_dl):
    row = cat_dl.infer('17:41197819-41199659:-', 'NA00002', 'psi5')
    assert row == {
        'junction': '17:41197819-41199659:-',
        'sample': 'NA00002',
        'count_cat': 12,
        'psi_cat': 1,
        'ref_psi_cat': 1.0,
        'k_cat': 65,
        'n_cat': 65,
        'delta_logit_psi_cat': 0,
        'delta_psi_cat': 0
    }

    row = cat_dl.infer('17:41251897-41256138:-', 'NA00001', 'psi5')
    assert row == {
        'junction': '17:41251897-41256138:-',
        'sample': 'NA00001',
        'count_cat': 9,
        'psi_cat': 0.6,
        'ref_psi_cat': 0.625,
        'k_cat': 30,
        'n_cat': 48,
        'delta_logit_psi_cat': -0.10536051565782634,
        'delta_psi_cat': -0.0205021934541626
    }
