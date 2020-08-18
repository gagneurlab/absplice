import pytest
from splicing_outlier_prediction import CatDataloader
from conftest import ref_table5_kn_file, ref_table3_kn_file, fasta_file, vcf_file, multi_vcf_file, count_cat_file


@pytest.fixture
def cat_dl3():
    return CatDataloader(
        ref_table3=ref_table3_kn_file,
        count_cat=count_cat_file)

@pytest.fixture
def cat_dl5():
    return CatDataloader(
        ref_table5=ref_table5_kn_file,
        count_cat=count_cat_file)


def test_cat_dataloader_init(cat_dl):
    assert cat_dl.ref_table5.method == 'kn'
    assert cat_dl.ref_table5.df.shape[0] == 26

    assert cat_dl.ref_table3.method == 'kn'
    assert cat_dl.ref_table3.df.shape[0] == 28


def test_cat_dataloader_iter(cat_dl):
    assert sum(1 for _ in cat_dl) > 0


def test_cat_dataloader_iter5(cat_dl5):
    assert sum(1 for _ in cat_dl5) > 0


def test_cat_dataloader_iter3(cat_dl3):
    assert sum(1 for _ in cat_dl3) > 0


def test_cat_dataloader_len(cat_dl):
    assert len(cat_dl) == 52
    
    
def test_cat_dataloader_len3(cat_dl3):
    assert len(cat_dl3) == 27
    
    
def test_cat_dataloader_len5(cat_dl5):
    assert len(cat_dl5) == 25
    
    
def test_cat_dataloader_getitem(cat_dl):
    rows = list(cat_dl) 
    rows = sorted(rows, key= lambda x: x['metadata']['target_tissue']['junction'])
    assert sorted(cat_dl.samples) == ['NA00001', 'NA00002', 'NA00003']
#     # OR
#     rows = [
#         i
#         for i in rows
#         if i['metadata']['target_tissue']['junction'] == '17:41201211-41203079:-'
#     ]

#     __import__("pdb").set_trace()
    
    assert rows[0]['metadata']['target_tissue']['junction'] == rows[0]['metadata']['cat_tissue']['junction'] 
    assert rows[0]['metadata']['target_tissue']['junction'] == '17:41197819-41199659:-' #first common junction, dataloader iterates over junctions
    assert rows[0]['metadata']['target_tissue']['event_type'] == 'psi5'
    assert rows[0]['metadata']['target_tissue']['psi'] == 1.0
    assert rows[0]['metadata']['cat_tissue']['counts'] == [25, 12, 28]
    assert rows[0]['metadata']['cat_tissue']['psi'] == [1.0, 1.0, 1.0]
    assert rows[0]['metadata']['cat_tissue']['psi_ref'] == 1.0
    assert rows[0]['metadata']['cat_tissue']['k'] == 65.0
    assert rows[0]['metadata']['cat_tissue']['n'] == 65.0
    
    assert rows[-1]['metadata']['target_tissue']['junction'] == rows[-1]['metadata']['cat_tissue']['junction'] 
    assert rows[-1]['metadata']['target_tissue']['junction'] == '17:41251897-41256138:-' #first common junction, dataloader iterates over junctions
    assert rows[-1]['metadata']['target_tissue']['event_type'] == 'psi5'
    assert rows[-1]['metadata']['target_tissue']['psi'] == 0.7454834226722256
    assert rows[-1]['metadata']['cat_tissue']['counts'] == [9, 6, 15]
    assert rows[-1]['metadata']['cat_tissue']['psi'] == [0.6, 0.6, 0.6521739130434783]
    assert rows[-1]['metadata']['cat_tissue']['psi_ref'] == 0.625
    assert rows[-1]['metadata']['cat_tissue']['k'] == 30.0
    assert rows[-1]['metadata']['cat_tissue']['n'] == 48.0
    
    
    