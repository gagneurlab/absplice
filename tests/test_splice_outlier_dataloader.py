import pytest
from splicing_outlier_prediction import SpliceOutlierDataloader
from conftest import ref_table5_kn_file, ref_table3_kn_file, ref_table5_kn_file2, ref_table3_kn_file2, fasta_file, vcf_file, multi_vcf_file, count_cat_file, combined_ref_tables5_file, combined_ref_tables3_file


@pytest.fixture
def outlier_dl5():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file, ref_tables5=[ref_table5_kn_file, ref_table5_kn_file2],
        combined_ref_tables5=combined_ref_tables5_file)


@pytest.fixture
def outlier_dl3():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file, ref_tables3=[ref_table3_kn_file])


def test_splicing_outlier_dataloader_init(outlier_dl):
    assert outlier_dl.combined_ref_tables5.method == ['kn', 'kn']
    assert outlier_dl.combined_ref_tables5.df.shape[0] == 50
    assert list(set(outlier_dl.combined_ref_tables5.ref_tables[0]['tissue'])) == ['lymphocytes']
    assert list(set(outlier_dl.combined_ref_tables5.ref_tables[1]['tissue'])) == ['lung']
    assert sorted(list(
        set(outlier_dl.combined_ref_tables5.ref_tables[0]['junctions']).union(
            set(outlier_dl.combined_ref_tables5.ref_tables[1]['junctions'])
        )
        )) == sorted(list(set(outlier_dl.combined_ref_tables5.df.index)))
    assert outlier_dl.combined_ref_tables3.method == ['kn', 'kn']
    assert outlier_dl.combined_ref_tables3.df.shape[0] == 51
    assert list(set(outlier_dl.combined_ref_tables3.ref_tables[0]['tissue'])) == ['lymphocytes']
    assert list(set(outlier_dl.combined_ref_tables3.ref_tables[1]['tissue'])) == ['lung']
    assert sorted(list(
        set(outlier_dl.combined_ref_tables3.ref_tables[0]['junctions']).union(
            set(outlier_dl.combined_ref_tables3.ref_tables[1]['junctions'])
        )
        )) == sorted(list(set(outlier_dl.combined_ref_tables3.df.index)))


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
    # assert junction['ref_psi'] == 1

    junction = rows[-1]['metadata']['junction']
    variant = rows[-1]['metadata']['variant']
    assert junction['junction'] == '17:41246877-41251791:-'
    assert variant['annotation'] == '17:41251886:A>G'
    assert junction['event_type'] == 'psi3'
    # assert junction['ref_psi'] == 0.2066666666666666
