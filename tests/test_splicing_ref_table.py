import pytest
from splicing_outlier_prediction import SplicingRefTable
from conftest import ref_table5_kn_file, ref_table3_kn_file


@pytest.fixture
def ref_table5_kn():
    return SplicingRefTable.read_csv([ref_table5_kn_file])#, regex_pattern='test_(.*)_ref')


@pytest.fixture
def ref_table3_kn():
    return SplicingRefTable.read_csv([ref_table3_kn_file], regex_pattern='test_(.*)_ref')


def test_SplicingRefTable__init__(ref_table5_kn, ref_table3_kn):
    assert ref_table5_kn.method == ['kn']
    assert ref_table5_kn.df.shape[0] == 26

    assert ref_table3_kn.method == ['kn']
    assert ref_table3_kn.df.shape[0] == 28


def test_SplicingRefTable_valid(ref_table5_kn, ref_table3_kn):
    pass
