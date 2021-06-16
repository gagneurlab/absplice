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
    assert outlier_dl.combined_ref_tables5.method == ['kn', 'kn']
    assert outlier_dl.combined_ref_tables5.df.shape[0] == 38
    assert list(set(outlier_dl.combined_ref_tables5.ref_tables[0]['tissue'])) == [
        'testis']
    assert list(set(outlier_dl.combined_ref_tables5.ref_tables[1]['tissue'])) == [
        'lung']
    assert sorted(list(
        set(outlier_dl.combined_ref_tables5.ref_tables[0]['junctions']).union(
            set(outlier_dl.combined_ref_tables5.ref_tables[1]['junctions'])
        )
    )) == sorted(list(set(outlier_dl.combined_ref_tables5.df.index)))
    assert outlier_dl.combined_ref_tables3.method == ['kn', 'kn']
    assert outlier_dl.combined_ref_tables3.df.shape[0] == 58
    assert list(set(outlier_dl.combined_ref_tables3.ref_tables[0]['tissue'])) == [
        'testis']
    assert list(set(outlier_dl.combined_ref_tables3.ref_tables[1]['tissue'])) == [
        'lung']
    assert sorted(list(
        set(outlier_dl.combined_ref_tables3.ref_tables[0]['junctions']).union(
            set(outlier_dl.combined_ref_tables3.ref_tables[1]['junctions'])
        )
    )) == sorted(list(set(outlier_dl.combined_ref_tables3.df.index)))


def test_splicing_outlier_dataloader_init_tissue_list():
    outlier_dl_2 = SpliceOutlierDataloader(
        fasta_file, vcf_file,
        splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
        splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
    )

    assert outlier_dl_2.combined_ref_tables5.method == ['kn', 'kn']
    assert outlier_dl_2.combined_ref_tables5.df.shape[0] == 38
    assert list(set(outlier_dl_2.combined_ref_tables5.ref_tables[0]['tissue'])) == [
        'testis']
    assert list(set(outlier_dl_2.combined_ref_tables5.ref_tables[1]['tissue'])) == [
        'lung']
    assert sorted(list(
        set(outlier_dl_2.combined_ref_tables5.ref_tables[0]['junctions']).union(
            set(outlier_dl_2.combined_ref_tables5.ref_tables[1]['junctions'])
        )
    )) == sorted(list(set(outlier_dl_2.combined_ref_tables5.df.index)))
    assert outlier_dl_2.combined_ref_tables3.method == ['kn', 'kn']
    assert outlier_dl_2.combined_ref_tables3.df.shape[0] == 58
    assert list(set(outlier_dl_2.combined_ref_tables3.ref_tables[0]['tissue'])) == [
        'testis']
    assert list(set(outlier_dl_2.combined_ref_tables3.ref_tables[1]['tissue'])) == [
        'lung']
    assert sorted(list(
        set(outlier_dl_2.combined_ref_tables3.ref_tables[0]['junctions']).union(
            set(outlier_dl_2.combined_ref_tables3.ref_tables[1]['junctions'])
        )
    )) == sorted(list(set(outlier_dl_2.combined_ref_tables3.df.index)))


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
    # assert junction['ref_psi'] == 1

    junction = rows[-1]['metadata']['junction']
    variant = rows[-1]['metadata']['variant']
    assert junction['junction'] == '17:41223255-41228504:-'
    assert variant['annotation'] == '17:41228514:CCTGGTTCTTTATTTTTACTGGT>C'
    assert junction['event_type'] == 'psi3'
    # assert junction['ref_psi'] == 0.2066666666666666
