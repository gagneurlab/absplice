import pytest
from splicing_outlier_prediction import SpliceOutlierDataloader
from conftest import ref_table5_kn_file, ref_table3_kn_file, fasta_file, vcf_file, multi_vcf_file, count_cat_file


@pytest.fixture
def outlier_dl5():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file, ref_table5=ref_table5_kn_file)


@pytest.fixture
def outlier_dl3():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file, ref_table3=ref_table3_kn_file)

@pytest.fixture
def outlier_dl_cat():
    return SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_table5=ref_table5_kn_file, ref_table3=ref_table3_kn_file, count_cat=count_cat_file,
        samples=True)


def test_splicing_outlier_dataloader_init(outlier_dl):
    assert outlier_dl.ref_table5.method == 'kn'
    assert outlier_dl.ref_table5.df.shape[0] == 26

    assert outlier_dl.ref_table3.method == 'kn'
    assert outlier_dl.ref_table3.df.shape[0] == 28


def test_splicing_outlier_dataloader_iter(outlier_dl):
    assert sum(1 for _ in outlier_dl) > 0


def test_splicing_outlier_dataloader_iter5(outlier_dl5):
    assert sum(1 for _ in outlier_dl5) > 0


def test_splicing_outlier_dataloader_iter3(outlier_dl3):
    assert sum(1 for _ in outlier_dl3) > 0

    
def test_splicing_outlier_dataloader_itercat(outlier_dl_cat):
    assert sum(1 for _ in outlier_dl3) > 0
    
    
def test_splicing_outlier_dataloader_next(outlier_dl):
    rows = list(outlier_dl)

    junction = rows[0]['metadata']['junction']
    variant = rows[0]['metadata']['variant']
    assert junction['junction'] == '17:41197819-41199659:-'
    assert variant['annotation'] == '17:41197805:ACATCTGCC>A'
    assert junction['event_type'] == 'psi5'
    assert junction['psi'] == 1

    junction = rows[-1]['metadata']['junction']
    variant = rows[-1]['metadata']['variant']
    assert junction['junction'] == '17:41246877-41251791:-'
    assert variant['annotation'] == '17:41251886:A>G'
    assert junction['event_type'] == 'psi3'
    assert junction['psi'] == 0.20666666666666667

    
def test_splicing_outlier_dataloader_cat_next(outlier_dl_cat):
    rows = list(outlier_dl_cat)

    junction = rows[0]['metadata']['junction']
    variant = rows[0]['metadata']['variant']
    
    assert junction['junction'] == '17:41201211-41203079:-' #for variant, take junction in 100bp window with highest mmsplice score
    assert variant['annotation'] == '17:41201201:TTC>CA' #from multi_vcf file
    assert junction['event_type'] == 'psi5'
    assert junction['psi'] == 0.9817811012916384

    junction = rows[-1]['metadata']['junction']
    variant = rows[-1]['metadata']['variant']
    assert junction['junction'] == '17:41199720-41201137:-'
    assert variant['annotation'] == '17:41201201:TTC>CA'
    assert junction['event_type'] == 'psi3'
    assert junction['psi'] == 0.9846668525229996
    
    cat = rows[0]['metadata']['junction']['cat']
    assert cat['samples'] == ['NA00002', 'NA00003']
    assert cat['counts'] == [12, 34]
    assert cat['psi'] == [1.0, 1.0]
    assert cat['psi_ref'] == 1.0
    assert cat['k'] == 83.0
    assert cat['n'] == 83.0
    