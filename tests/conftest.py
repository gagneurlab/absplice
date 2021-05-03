import pytest
from splicing_outlier_prediction import SpliceOutlierDataloader, SpliceOutlier, CatInference

vcf_file = 'tests/data/test.vcf.gz'
multi_vcf_file = 'tests/data/multi_test.vcf.gz'
fasta_file = 'tests/data/hg19.nochr.chr17.fa'
ref_table5_kn_testis = 'tests/data/test_testis_ref_table5_kn.csv'
ref_table3_kn_testis = 'tests/data/test_testis_ref_table3_kn.csv'
ref_table5_kn_lung = 'tests/data/test_lung_ref_table5_kn.csv'
ref_table3_kn_lung = 'tests/data/test_lung_ref_table3_kn.csv'
combined_ref_tables5_testis_lung = 'tests/data/test_combined_ref_tables5_testis_lung_kn.csv'
combined_ref_tables3_testis_lung = 'tests/data/test_combined_ref_tables3_testis_lung_kn.csv'
count_cat_file_lymphocytes = 'tests/data/test_count_table_cat_chrom17_lymphocytes.csv'
count_cat_file_blood = 'tests/data/test_count_table_cat_chrom17_blood.csv'
spliceAI = 'tests/data/test_spliceAI.csv'
pickle_DNA = 'tests/data/model_DNA_trained_on_all_GTEx.pkl'
pickle_DNA_CAT = 'tests/data/model_CAT_concat.pkl'


@pytest.fixture
def outlier_dl():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file,
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung], 
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        combined_ref_tables5=combined_ref_tables5_testis_lung, 
        combined_ref_tables3=combined_ref_tables3_testis_lung,
        regex_pattern='test_(.*)_ref',
        save_combined_ref_tables=True,
        )

@pytest.fixture
def cat_dl():
    return CatInference(ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung], 
                        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
                        regex_pattern='test_(.*)_ref',
                        count_cat=[count_cat_file_lymphocytes, count_cat_file_blood],
                        regex_pattern_cat='chrom17_(.*).csv',
                        )

@pytest.fixture
def outlier_model():
    return SpliceOutlier()

@pytest.fixture
def outlier_results(outlier_model, outlier_dl):
    return outlier_model.predict_on_dataloader(outlier_dl)


@pytest.fixture
def outlier_dl_multi():
    return SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_testis, ref_table5_kn_lung], 
        ref_tables3=[ref_table3_kn_testis, ref_table3_kn_lung],
        combined_ref_tables5=combined_ref_tables5_testis_lung, 
        combined_ref_tables3=combined_ref_tables3_testis_lung,
        regex_pattern='test_(.*)_ref',
        samples=True
        )
