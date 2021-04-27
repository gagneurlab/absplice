import pytest
from splicing_outlier_prediction import SplicingRefTable
from conftest import ref_table5_kn_testis, ref_table3_kn_testis


@pytest.fixture
def ref_table5_kn():
    return SplicingRefTable.read_csv([ref_table5_kn_testis])#, regex_pattern='test_(.*)_ref')


@pytest.fixture
def ref_table3_kn():
    return SplicingRefTable.read_csv([ref_table3_kn_testis], regex_pattern='test_(.*)_ref')


def test_SplicingRefTable__init__(ref_table5_kn, ref_table3_kn):
    assert ref_table5_kn.method == ['kn']
    assert ref_table5_kn.df.shape[0] == 24

    assert ref_table3_kn.method == ['kn']
    assert ref_table3_kn.df.shape[0] == 45


def test_SplicingRefTable_valid(ref_table5_kn, ref_table3_kn):
    pass
