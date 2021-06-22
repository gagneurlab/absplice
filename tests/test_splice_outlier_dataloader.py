import pytest
from splicing_outlier_prediction import SpliceOutlierDataloader
from conftest import ref_table5_kn_testis, ref_table3_kn_testis, ref_table5_kn_lung, ref_table3_kn_lung, fasta_file, vcf_file, multi_vcf_file, count_cat_file_lymphocytes, combined_ref_tables5_testis_lung, combined_ref_tables3_testis_lung


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
    assert outlier_dl.combined_splicemap5.shape[0] == 38

    assert outlier_dl.splicemaps5[0].name == 'gtex-grch37-testis'
    assert outlier_dl.splicemaps5[1].name == 'gtex-grch37-lung'

    assert sorted(
        set(outlier_dl.splicemaps5[0].df['junctions']).union(
            set(outlier_dl.splicemaps5[1].df['junctions'])
        )
    ) == sorted(set(outlier_dl.combined_splicemap5.index))

    # splicemap3
    assert outlier_dl.combined_splicemap3.shape[0] == 58

    assert outlier_dl.splicemaps3[0].name == 'gtex-grch37-testis'
    assert outlier_dl.splicemaps3[1].name == 'gtex-grch37-lung'

    assert sorted(
        set(outlier_dl.splicemaps3[0].df['junctions']).union(
            set(outlier_dl.splicemaps3[1].df['junctions'])
        )
    ) == sorted(set(outlier_dl.combined_splicemap3.index))


def test_splicing_outlier_dataloader_init_dl3(outlier_dl3):
    # psi 5
    assert outlier_dl3.combined_splicemap5 is None
    # psi 3
    assert outlier_dl3.combined_splicemap3.shape[0] == 45
    assert outlier_dl3.splicemaps3[0].name == 'gtex-grch37-testis'

    assert sorted(set(outlier_dl3.splicemaps3[0].df['junctions'])) \
        == sorted(set(outlier_dl3.combined_splicemap3.index))


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
    assert junction['junction'] == '17:41199720-41201137:-'
    assert variant['annotation'] == '17:41199651:G>A'
    assert junction['event_type'] == 'psi5'

    junction = rows[-1]['metadata']['junction']
    variant = rows[-1]['metadata']['variant']
    assert junction['junction'] == '17:41223255-41228504:-'
    assert variant['annotation'] == '17:41228514:CCTGGTTCTTTATTTTTACTGGT>C'
    assert junction['event_type'] == 'psi3'
