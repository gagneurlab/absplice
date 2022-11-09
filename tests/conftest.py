import pytest
import pandas as pd
import tempfile
from absplice import SplicingOutlierResult, \
    SpliceOutlierDataloader, SpliceOutlier, CatInference, \
        GENE_MAP, GENE_TPM

vcf_file = 'tests/data/test.vcf.gz'
multi_vcf_file = 'tests/data/multi_test.vcf.gz'
var_samples_path = 'tests/data/multi_test.vcf_samples.csv'
fasta_file = 'tests/data/hg19.nochr.chr17.fa'

ref_table5_kn_testis = 'tests/data/Testis_splicemap_psi5_method=kn_event_filter=median_cutoff.csv.gz'
ref_table3_kn_testis = 'tests/data/Testis_splicemap_psi3_method=kn_event_filter=median_cutoff.csv.gz'
ref_table5_kn_lung = 'tests/data/Lung_splicemap_psi5_method=kn_event_filter=median_cutoff.csv.gz'
ref_table3_kn_lung = 'tests/data/Lung_splicemap_psi3_method=kn_event_filter=median_cutoff.csv.gz'
ref_table5_kn_blood = 'tests/data/Cells_Cultured_fibroblasts_splicemap_psi5_method=kn_event_filter=median_cutoff.csv.gz'
ref_table3_kn_blood = 'tests/data/Cells_Cultured_fibroblasts_splicemap_psi3_method=kn_event_filter=median_cutoff.csv.gz'

count_cat_file_lymphocytes_complete = 'tests/data/create_test/data/full_data/backup/test_count_table_cat_chrom17_lymphocytes.csv'
count_cat_file_blood_complete = 'tests/data/create_test/data/full_data/backup/test_count_table_cat_chrom17_blood.csv'
count_cat_file_lymphocytes = 'tests/data/test_count_table_cat_chrom17_lymphocytes.csv'
count_cat_file_blood = 'tests/data/test_count_table_cat_chrom17_blood.csv'

mmsplice_path = 'tests/data/test_mmsplice.csv'
spliceai_path = 'tests/data/test_spliceAI.csv'
mmsplice_cat_path = 'tests/data/test_mmsplice_cat.csv'

spliceai_vcf_path = 'tests/data/spliceai_snv.vcf'
spliceai_vcf_path2 = 'tests/data/test_spliceai.vcf'

cadd_splice_path = 'tests/data/cadd_splice_test.tsv.gz'

absplice_precomputed_path = 'tests/data/all_SNVs_precomputed.csv'

@pytest.fixture
def df_var_samples():
    return pd.read_csv(var_samples_path)

@pytest.fixture
def df_mmsplice():
    return pd.read_csv(mmsplice_path)

@pytest.fixture
def df_spliceai():
    return pd.read_csv(spliceai_path)

@pytest.fixture
def df_mmsplice_cat():
    return pd.read_csv(mmsplice_cat_path)

@pytest.fixture
def outlier_dl():
    return SpliceOutlierDataloader(
        fasta_file, vcf_file,
        splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
        splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung])

@pytest.fixture
def outlier_dl_multi():
    return SpliceOutlierDataloader(
        fasta_file, multi_vcf_file,
        splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
        splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung]
    )

@pytest.fixture
def cat_dl():
    return [
        CatInference(splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
                     splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
                     count_cat=count_cat_file_lymphocytes,
                     name='lymphocytes'),
        CatInference(splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
                     splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
                     count_cat=count_cat_file_blood,
                     name='blood'),
    ]

@pytest.fixture
def outlier_model():
    return SpliceOutlier()

@pytest.fixture
def outlier_results(outlier_model, outlier_dl):
    return outlier_model.predict_on_dataloader(outlier_dl)

@pytest.fixture
def outlier_results_multi(outlier_model, outlier_dl_multi, df_var_samples):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    results.add_samples(df_var_samples)
    return results

@pytest.fixture
def gene_map():
    return pd.read_csv(GENE_MAP, sep='\t')

@pytest.fixture
def gene_tpm():
    return pd.read_csv(GENE_TPM)

# @pytest.fixture
# def outlier_results_complete(outlier_model, outlier_dl_multi, df_spliceai, gene_map, df_var_samples):
#     results = outlier_model.predict_on_dataloader(outlier_dl_multi)
#     results.add_spliceai(df_spliceai, gene_map)
#     results.add_samples(df_var_samples)
#     return results

@pytest.fixture
def mmsplice_splicemap_cols():
    return sorted([
            'variant', 'tissue', 'junction', 'event_type',
            'splice_site', 'ref_psi', 'median_n', 
            'gene_id', 'gene_name', 'gene_tpm',
            'delta_logit_psi', 'delta_psi',
        ])
     
variants = [
    "chr3:193360794:C:['A']"
]

def parse_vcf_id(vcf_id):
    return vcf_id.replace("'", '').replace('[', '').replace(']', '').split(':')

@pytest.fixture
def vcf_path():
    with tempfile.NamedTemporaryFile('w') as temp_vcf:
        temp_vcf.write('##fileformat=VCFv4.0\n')
        temp_vcf.write('##contig=<ID=13,length=115169878>\n')
        temp_vcf.write('##contig=<ID=17,length=81195210>\n')
        temp_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        for v in variants:
            temp_vcf.write('%s\t%s\t1\t%s\t%s\t.\t.\t.\n'
                           % tuple(parse_vcf_id(v)))

        temp_vcf.flush()
        yield temp_vcf.name