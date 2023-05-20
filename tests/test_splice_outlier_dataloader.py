import pytest
from absplice import SpliceOutlierDataloader
from conftest import fasta_file, vcf_file, multi_vcf_file, \
    ref_table5_kn_testis, ref_table3_kn_testis, ref_table5_kn_lung, ref_table3_kn_lung, \
    count_cat_file_lymphocytes

from absplice.utils import junction_str_to_tuple


@pytest.fixture
def outlier_dl5():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file,
        splicemap5=[
            ref_table5_kn_testis,
            ref_table5_kn_lung
        ])


@pytest.fixture
def outlier_dl3():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file, splicemap3=ref_table3_kn_testis)


def test_splicing_outlier_dataloader_init(outlier_dl):
    # splicemap5
    assert outlier_dl.combined_splicemap5.shape[0] == 6

    assert outlier_dl.splicemaps5[0].name == 'Testis'
    assert outlier_dl.splicemaps5[1].name == 'Lung'

    assert [junction_str_to_tuple(k) for k in sorted(
        set(outlier_dl.splicemaps5[0].df['junctions']).union(
            set(outlier_dl.splicemaps5[1].df['junctions'])
        )
    )] == list(sorted(set(outlier_dl.combined_splicemap5)))

    # splicemap3
    assert outlier_dl.combined_splicemap3.shape[0] == 5

    assert outlier_dl.splicemaps3[0].name == 'Testis'
    assert outlier_dl.splicemaps3[1].name == 'Lung'

    assert [junction_str_to_tuple(k) for k in sorted(
        set(outlier_dl.splicemaps3[0].df['junctions']).union(
            set(outlier_dl.splicemaps3[1].df['junctions'])
        )
    )] == list(sorted(set(outlier_dl.combined_splicemap3)))


def test_splicing_outlier_dataloader_feather(outlier_dl5_feather):
    # psi 5
    assert outlier_dl5_feather.combined_splicemap5 is not None
    assert len(outlier_dl5_feather.combined_splicemap5) == 5
    assert len(outlier_dl5_feather.splicemaps5) == 1
    assert outlier_dl5_feather.splicemaps5[0].name == 'Testis'
    # psi 3
    assert outlier_dl5_feather.combined_splicemap3 is None


def test_splicing_outlier_dataloader_init_dl3(outlier_dl3):
    # psi 5
    assert outlier_dl3.combined_splicemap5 is None
    # psi 3
    assert outlier_dl3.combined_splicemap3.shape[0] == 4
    assert outlier_dl3.splicemaps3[0].name == 'Testis'

    assert (
            [junction_str_to_tuple(k) for k in sorted(set(outlier_dl3.splicemaps3[0].df['junctions']))]
            == list(sorted(set(outlier_dl3.combined_splicemap3)))
    )


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
    # assert junction['junction'] == '17:41201211-41203079:-'
    assert junction['junction'] == '17:41201200-41205000:-'
    assert variant['annotation'] == '17:41201201:TTC>CA'
    assert junction['event_type'] == 'psi5'

    junction = rows[-1]['metadata']['junction']
    variant = rows[-1]['metadata']['variant']
    assert junction['junction'] == '17:41267796-41276033:-'
    assert variant['annotation'] == '17:41276032:T>A'
    assert junction['event_type'] == 'psi3'
