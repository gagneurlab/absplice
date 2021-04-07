import pytest
from splicing_outlier_prediction import SpliceOutlierDataloader
from conftest import ref_table5_kn_file, ref_table3_kn_file, ref_table5_kn_file2, ref_table3_kn_file2, fasta_file, vcf_file, multi_vcf_file, count_cat_file, intron_annotation5_file, intron_annotation3_file


@pytest.fixture
def outlier_dl5():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file, ref_tables5=[ref_table5_kn_file, ref_table5_kn_file2],
        intron_annotation5=intron_annotation5_file)


@pytest.fixture
def outlier_dl3():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file, ref_tables3=[ref_table3_kn_file])


def test_splicing_outlier_dataloader_init(outlier_dl):
    assert outlier_dl.intron_annotation5.method == 'kn'
    assert outlier_dl.intron_annotation5.df.shape[0] == 50

    assert outlier_dl.intron_annotation3.method == 'kn'
    assert outlier_dl.intron_annotation3.df.shape[0] == 51


def test_splicing_outlier_dataloader_iter(outlier_dl):
    assert sum(1 for _ in outlier_dl) > 0


def test_splicing_outlier_dataloader_iter5(outlier_dl5):
    assert sum(1 for _ in outlier_dl5) > 0


def test_splicing_outlier_dataloader_iter3(outlier_dl3):
    assert sum(1 for _ in outlier_dl3) > 0


def test_splicing_outlier_dataloader_next(outlier_dl):
    rows = list(outlier_dl)

    junction = rows[0]['metadata']['junction']
    variant = rows[0]['metadata']['variant']
    assert junction['junction'] == '17:41197819-41199659:-'
    assert variant['annotation'] == '17:41197805:ACATCTGCC>A'
    assert junction['event_type'] == 'psi5'
    assert junction['ref_psi'] == 1

    junction = rows[-1]['metadata']['junction']
    variant = rows[-1]['metadata']['variant']
    assert junction['junction'] == '17:41246877-41251791:-'
    assert variant['annotation'] == '17:41251886:A>G'
    assert junction['event_type'] == 'psi3'
    assert junction['ref_psi'] == 0.2066666666666666
