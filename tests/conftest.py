import pytest
from splicing_outlier_prediction import SpliceOutlierDataloader


vcf_file = 'tests/data/test.vcf.gz'
multi_vcf_file = 'tests/data/multi_test.vcf.gz'
fasta_file = 'tests/data/hg19.nochr.chr17.fa'
ref_table5_kn_file = 'tests/data/test_lymphocytes_ref_table5_kn.csv'
ref_table3_kn_file = 'tests/data/test_lymphocytes_ref_table3_kn.csv'
count_cat_file = 'tests/data/test_count_table_cat_chrom17.csv'


@pytest.fixture
def outlier_dl():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file,
        ref_table5=ref_table5_kn_file, ref_table3=ref_table3_kn_file,
        count_cat=count_cat_file)
