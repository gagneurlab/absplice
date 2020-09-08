import pytest
import numpy as np
from splicing_outlier_prediction import CatDataloader
from conftest import ref_table5_kn_file, ref_table3_kn_file, fasta_file, \
    vcf_file, multi_vcf_file, count_cat_file


@pytest.fixture
def cat_dl5():
    return CatDataloader(
        ref_table5=ref_table5_kn_file,
        count_cat=count_cat_file)


@pytest.fixture
def cat_dl3():
    return CatDataloader(
        ref_table3=ref_table3_kn_file,
        count_cat=count_cat_file)


def test_cat_dataloader_init(cat_dl):
    assert cat_dl.ref_table5.method == 'kn'
    assert cat_dl.ref_table5.df.shape[0] == 26

    assert cat_dl.ref_table3.method == 'kn'
    assert cat_dl.ref_table3.df.shape[0] == 28


def test_cat_dataloader_iter(cat_dl):
    assert sum(1 for _ in cat_dl) == 52


def test_cat_dataloader_iter5(cat_dl5):
    assert sum(1 for _ in cat_dl5) == 25


def test_cat_dataloader_iter3(cat_dl3):
    assert sum(1 for _ in cat_dl3) == 27


def test_cat_dataloader_len(cat_dl):
    assert len(cat_dl) == 52


def test_cat_dataloader_len3(cat_dl3):
    assert len(cat_dl3) == 27


def test_cat_dataloader_len5(cat_dl5):
    assert len(cat_dl5) == 25


def test_cat_dataloader_getitem(cat_dl):
    rows = list(cat_dl)
    rows = sorted(rows, key=lambda x: x['metadata']['junction'])
    assert sorted(cat_dl.samples) == ['NA00001', 'NA00002', 'NA00003']

    np.testing.assert_almost_equal(
        rows[0]['inputs']['counts'],
        [25, 12, 28]
    )
    np.testing.assert_almost_equal(
        rows[0]['inputs']['psi'],
        [1., 1., 1.]
    )
    assert rows[0]['metadata'] == {
        'junction': '17:41197819-41199659:-',
        'event_type': 'psi5',
        'cat_tissue': {
            'ref_psi': 1.0,
            'k': 65,
            'n': 65
        },
        'target_tissue': {
            'Chromosome': 17,
            'Start': 41197819,
            'End': 41199659,
            'Strand': '-',
            'events': '17:41197819-41199659:-',
            'splice_site': '17:41199659:-',
            'ref_psi': 1.0, 'k': 6337, 'n': 6337,
            'gene_id': 'ENSG00000012048.22_5', 'gene_name': 'BRCA1',
            'weak': False,
            'transcript_id': 'ENST00000586385.5_1;ENST00000591534.5_1;'
            'ENST00000461221.5_1;ENST00000493795.5_1;ENST00000357654.8_3;'
            'ENST00000591849.5_1;ENST00000468300.5_2;ENST00000471181.7_3;'
            'ENST00000491747.6_3;ENST00000352993.7_2;ENST00000644379.1_1',
            'gene_type': 'protein_coding'
        }
    }

    np.testing.assert_almost_equal(
        rows[-1]['inputs']['counts'],
        [9, 6, 15]
    )
    np.testing.assert_almost_equal(
        rows[-1]['inputs']['psi'],
        [0.6, 0.6, 0.6521739130434783]
    )
    assert rows[-1]['metadata'] == {
        'junction': '17:41251897-41256138:-',
        'event_type': 'psi5',
        'cat_tissue': {
            'ref_psi': 0.625,
            'k': 30,
            'n': 48
        },
        'target_tissue': {
            'Chromosome': 17,
            'Start': 41251897,
            'End': 41256138,
            'Strand': '-',
            'events': '17:41251894-41256138:-;17:41251897-41256138:-',
            'splice_site': '17:41256138:-',
            'ref_psi': 0.7454834226722256, 'k': 3755, 'n': 5037,
            'gene_id': 'ENSG00000012048.22_5', 'gene_name': 'BRCA1',
            'weak': False,
            'transcript_id': 'ENST00000470026.5_1;ENST00000493919.5_1;'
            'ENST00000354071.7_1;ENST00000634433.1_1;ENST00000494123.5_1;'
            'ENST00000493795.5_1;ENST00000357654.8_3;ENST00000487825.5_1;'
            'ENST00000468300.5_2;ENST00000471181.7_3;ENST00000477152.5_2;'
            'ENST00000492859.5_1;ENST00000642945.1_1;ENST00000652672.1_1;'
            'ENST00000491747.6_3;ENST00000352993.7_2;ENST00000461798.5_1',
            'gene_type': 'protein_coding'
        }
    }
