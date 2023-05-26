import pytest
import pandas as pd
import numpy as np
from absplice import SpliceOutlier, SpliceOutlierDataloader, CatInference
from absplice.ensemble import train_model_ebm
from conftest import fasta_file, multi_vcf_file, \
    ref_table5_kn_testis, ref_table3_kn_testis,  \
    ref_table5_kn_lung, ref_table3_kn_lung, \
    count_cat_file_lymphocytes,  count_cat_file_blood, \
    spliceai_path


def test_splicing_outlier_on_batch(outlier_model, outlier_dl, mmsplice_splicemap_cols):
    batch = next(outlier_dl.batch_iter())
    df = outlier_model.predict_on_batch(batch, outlier_dl)
    assert sorted(df.columns.tolist()) == mmsplice_splicemap_cols
    assert df.shape[0] == 4


def test_splicing_outlier_predict_on_dataloader(outlier_model, outlier_dl, mmsplice_splicemap_cols):
    results = outlier_model.predict_on_dataloader(outlier_dl)
    assert sorted(results.df_mmsplice.columns.tolist()) == mmsplice_splicemap_cols


def test_splicing_outlier_predict_save(outlier_model, outlier_dl, tmp_path, mmsplice_splicemap_cols):
    # output_csv = '/home/wagnern/Projects/pred.csv'
    output_csv = tmp_path / 'test_mmsplice.csv'
    outlier_model.predict_save(outlier_dl, output_csv)
    df = pd.read_csv(output_csv)
    assert sorted(df.columns.tolist()) == mmsplice_splicemap_cols


def test_multi_sample_predict(outlier_dl_multi, outlier_model):
    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    print(results.df_mmsplice.shape)
    
    
def test_var_negative_strand(vcf_path):
    from conftest import ref_table3_dummy, ref_table3_dummyTESTIS, ref_table5_dummy
    dl = SpliceOutlierDataloader(
        fasta_file, vcf_path,
        splicemap5=[ref_table5_dummy],
        splicemap3=[ref_table3_dummyTESTIS, ref_table3_dummy]
    )
    model = SpliceOutlier()
    
    results = model.predict_on_dataloader(dl)
    df = results.df_mmsplice
    assert df[df['variant'].str.contains('17:41203228:T>A')]['delta_logit_psi'].max() < 10
    
    
def test_singleVariantMatcher(vcf_path):
    from kipoiseq.extractors import MultiSampleVCF, SingleVariantMatcher
    import pyranges as pr
    pr_exons = pd.read_csv('tests/data/pr_exons_debug.csv')
    pr_exons = pr.PyRanges(pr_exons)
    
    matcher = SingleVariantMatcher(
        vcf_path, pranges=pr_exons,
        interval_attrs=('junction',)
    )
    
    rows = list(matcher)
    assert len(rows) == 2
    
    
    
