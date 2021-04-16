import pytest
from splicing_outlier_prediction import SpliceOutlierDataloader, CatInference


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


@pytest.fixture
def outlier_dl():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file,
        ref_tables5=[ref_table5_kn_file, ref_table5_kn_file2], ref_tables3=[ref_table3_kn_file, ref_table3_kn_file2],
        combined_ref_tables5=combined_ref_tables5_file, 
        combined_ref_tables3=combined_ref_tables3_file,
        regex_pattern='test_(.*)_ref'
        )


@pytest.fixture
def cat_dl():
    return CatInference(ref_tables5=[ref_table5_kn_file],
                        ref_tables3=[ref_table3_kn_file],
                        count_cat=count_cat_file)
