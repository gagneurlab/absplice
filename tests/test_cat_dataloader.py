import pytest
from splicing_outlier_prediction import CatInference
from conftest import ref_table5_kn_testis, ref_table3_kn_testis, \
    ref_table5_kn_lung, ref_table3_kn_lung, \
    count_cat_file_lymphocytes, count_cat_file_blood, \
    vcf_file, multi_vcf_file, fasta_file


@pytest.fixture
def cat_dl5():
    return CatInference(
        ref_tables5=[ref_table5_kn_testis],
        regex_pattern='test_(.*)_ref',
        count_cat=[count_cat_file_lymphocytes],
        regex_pattern_cat='chrom17_(.*).csv'
    )


@pytest.fixture
def cat_dl3():
    return CatInference(
        ref_tables3=[ref_table3_kn_testis],
        regex_pattern='test_(.*)_ref',
        count_cat=[count_cat_file_lymphocytes],
        regex_pattern_cat='chrom17_(.*).csv'
    )


def test_cat_dataloader_init(cat_dl):
    assert cat_dl.combined_ref_tables5.method == ['kn', 'kn']
    assert cat_dl.combined_ref_tables5.df_all.shape[0] == 49
    assert sorted(list(cat_dl.combined_ref_tables5.df_all['tissue'].unique())) \
        == sorted(['testis', 'lung'])

    assert cat_dl.combined_ref_tables3.method == ['kn', 'kn']
    assert cat_dl.combined_ref_tables3.df_all.shape[0] == 71
    assert sorted(list(cat_dl.combined_ref_tables3.df_all['tissue'].unique())) \
        == sorted(['testis', 'lung'])


def test_cat_dataloader_common(cat_dl):
    assert len(cat_dl.common_junctions5[0]) == 35
    assert len(cat_dl.common_junctions3[0]) == 41


def test_cat_dataloader_common5(cat_dl5):
    assert len(cat_dl5.common_junctions5[0]) == 22


def test_cat_dataloader_common3(cat_dl3):
    assert len(cat_dl3.common_junctions3[0]) == 29


def test_cat_dataloader_sample_mapping():
    sample_mapping = {
        'NA00001': 'new_name1',
        'NA00002': 'new_name2',
        'NA00003': 'new_name3'
    }

    cat_dl_map = CatInference(
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung],
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        regex_pattern='test_(.*)_ref',
        count_cat=[count_cat_file_lymphocytes, count_cat_file_blood],
        regex_pattern_cat='chrom17_(.*).csv',
        sample_mapping=sample_mapping
    )

    assert sorted(list(cat_dl_map.samples[0])) \
        == sorted(list({'new_name1', 'new_name2', 'new_name3'}))


def test_cat_dataloader_infer(cat_dl):
    row = cat_dl.infer('17:41154800-41154888:+', 'NA00002', 'psi5')
    assert row == [
        {
            'junction': '17:41154800-41154888:+',
            'sample': 'NA00002',
            'tissue': 'testis',
            'tissue_cat': 'lymphocytes',
            'count_cat': 3936,
            'delta_logit_psi_cat': 0.0,
            'delta_psi_cat': 0.0,
            'k_cat': 15772,
            'median_n_cat': 5596.0,
            'n_cat': 15772,
            'psi_cat': 1.0,
            'ref_psi_cat': 1.0
        }, {
            'junction': '17:41154800-41154888:+',
            'sample': 'NA00002',
            'tissue': 'lung',
            'tissue_cat': 'lymphocytes',
            'count_cat': 3936,
            'delta_logit_psi_cat': 0.0,
            'delta_psi_cat': 0.0,
            'k_cat': 15772,
            'median_n_cat': 5596.0,
            'n_cat': 15772,
            'psi_cat': 1.0,
            'ref_psi_cat': 1.0,
        }, {
            'junction': '17:41154800-41154888:+',
            'sample': 'NA00002',
            'tissue': 'testis',
            'tissue_cat': 'blood',
            'count_cat': 438,
            'delta_logit_psi_cat': 0.0,
            'delta_psi_cat': 0.0,
            'k_cat': 1820,
            'median_n_cat': 502.0,
            'n_cat': 1820,
            'psi_cat': 1.0,
            'ref_psi_cat': 1.0,
        }, {
            'junction': '17:41154800-41154888:+',
            'sample': 'NA00002',
            'tissue': 'lung',
            'tissue_cat': 'blood',
            'count_cat': 438,
            'delta_logit_psi_cat': 0.0,
            'delta_psi_cat': 0.0,
            'k_cat': 1820,
            'median_n_cat': 502.0,
            'n_cat': 1820,
            'psi_cat': 1.0,
            'ref_psi_cat': 1.0,
        }]


def test_cat_dataloader_contains(cat_dl):
    assert cat_dl.contains('17:41201917-41203079:-', 'NA00002', 'psi5')


def test_cat_dataloader_contains_not(cat_dl):
    assert not cat_dl.contains(
        '17:41201917-4120307900000:-', 'NA00002', 'psi5')


def test_cat_dataloader_infer_only_one_cat(cat_dl):
    row = cat_dl.infer('17:41201917-41203079:-', 'NA00002', 'psi5')
    assert row == [{
        'junction': '17:41201917-41203079:-',
        'sample': 'NA00002',
        'tissue': 'testis',
        'tissue_cat': 'lymphocytes',
        'count_cat': 0,
        'delta_logit_psi_cat': 0.0,
        'delta_psi_cat': 1.734723475976807e-18,
        'k_cat': 0,
        'median_n_cat': 43.0,
        'n_cat': 166,
        'psi_cat': 0.0,
        'ref_psi_cat': 0.0
    }]
