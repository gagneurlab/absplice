import pytest
from splicing_outlier_prediction import SpliceOutlierDataloader, SpliceOutlier, CatInference

vcf_file = 'tests/data/test.vcf.gz'
multi_vcf_file = 'tests/data/multi_test.vcf.gz'
fasta_file = 'tests/data/hg19.nochr.chr17.fa'
ref_table5_kn_file = 'tests/data/test_lymphocytes_ref_table5_kn.csv'
ref_table3_kn_file = 'tests/data/test_lymphocytes_ref_table3_kn.csv'
ref_table5_kn_file2 = 'tests/data/test_lung_ref_table5_kn.csv'
ref_table3_kn_file2 = 'tests/data/test_lung_ref_table3_kn.csv'
combined_ref_tables5_file = 'tests/data/test_combined_ref_tables5_kn.csv'
combined_ref_tables3_file = 'tests/data/test_combined_ref_tables3_kn.csv'
count_cat_file = 'tests/data/test_count_table_cat_chrom17.csv'
spliceai_db_path = 'tests/data/spliceAI.db'

ref_table5_kn_file_fake1 = 'tests/data/test_brain_ref_table5_kn_fake.csv'
ref_table3_kn_file_fake1 = 'tests/data/test_brain_ref_table3_kn_fake.csv'
ref_table5_kn_file_fake2 = 'tests/data/test_liver_ref_table5_kn_fake.csv'
ref_table3_kn_file_fake2 = 'tests/data/test_liver_ref_table3_kn_fake.csv'
combined_ref_tables5_file_fake = 'tests/data/test_combined_ref_tables5_kn_fake.csv'
combined_ref_tables3_file_fake = 'tests/data/test_combined_ref_tables3_kn_fake.csv'


@pytest.fixture
def outlier_dl():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file,
        ref_tables5=[ref_table5_kn_file, ref_table5_kn_file2], 
        ref_tables3=[ref_table3_kn_file, ref_table3_kn_file2],
        combined_ref_tables5=combined_ref_tables5_file, 
        combined_ref_tables3=combined_ref_tables3_file,
        regex_pattern='test_(.*)_ref'
        )

@pytest.fixture
def cat_dl():
    return CatInference(ref_tables5=[ref_table5_kn_file, ref_table5_kn_file2], 
                        ref_tables3=[ref_table3_kn_file, ref_table3_kn_file2],
                        regex_pattern='test_(.*)_ref',
                        count_cat=count_cat_file)

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
        ref_tables5=[ref_table5_kn_file, ref_table5_kn_file2], 
        ref_tables3=[ref_table3_kn_file, ref_table3_kn_file2],
        combined_ref_tables5=combined_ref_tables5_file, 
        combined_ref_tables3=combined_ref_tables3_file,
        regex_pattern='test_(.*)_ref',
        samples=True
        )


@pytest.fixture
def outlier_dl_fake():
    return SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        ref_tables5=[ref_table5_kn_file_fake1, ref_table5_kn_file_fake2], 
        ref_tables3=[ref_table3_kn_file_fake1, ref_table3_kn_file_fake2],
        combined_ref_tables5=combined_ref_tables5_file_fake, 
        combined_ref_tables3=combined_ref_tables3_file_fake,
        regex_pattern='test_(.*)_ref',
        samples=True
        )

@pytest.fixture
def cat_dl_fake():
    return CatInference(ref_tables5=[ref_table5_kn_file_fake1, ref_table5_kn_file_fake2], 
                        ref_tables3=[ref_table3_kn_file_fake1, ref_table3_kn_file_fake2],
                        regex_pattern='test_(.*)_ref',
                        count_cat=count_cat_file)
