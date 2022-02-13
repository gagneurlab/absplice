import pytest
import pandas as pd
import numpy as np
from splicing_outlier_prediction import SpliceOutlier, SpliceOutlierDataloader, CatInference
from splicing_outlier_prediction.ensemble import train_model_ebm
from conftest import fasta_file, multi_vcf_file, \
    ref_table5_kn_testis, ref_table3_kn_testis,  \
    ref_table5_kn_lung, ref_table3_kn_lung, \
    count_cat_file_lymphocytes,  count_cat_file_blood, \
    spliceai_path, pickle_DNA, pickle_DNA_CAT, mmsplice_splicemap_cols


def test_splicing_outlier_on_batch(outlier_model, outlier_dl, mmsplice_splicemap_cols):
    batch = next(outlier_dl.batch_iter())
    df = outlier_model.predict_on_batch(batch, outlier_dl)
    assert sorted(df.columns.tolist()) == mmsplice_splicemap_cols
    assert df.shape[0] == 4


def test_splicing_outlier_predict_on_dataloader(outlier_model, outlier_dl, mmsplice_splicemap_cols):
    results = outlier_model.predict_on_dataloader(outlier_dl)
    assert sorted(results.df_mmsplice.columns.tolist()) == mmsplice_splicemap_cols


# def test_splicing_outlier_predict_on_dataloader_correct_tissue(outlier_model, outlier_dl):
#     results = outlier_model.predict_on_dataloader(outlier_dl)
#     assert False not in results.df.apply(lambda x: x['event_type'] in x['tissue'], axis=1)


def test_splicing_outlier_predict_save(outlier_model, outlier_dl, tmp_path, mmsplice_splicemap_cols):
    # output_csv = '/home/wagnern/Projects/pred.csv'
    output_csv = tmp_path / 'test_mmsplice.csv'
    outlier_model.predict_save(outlier_dl, output_csv)
    df = pd.read_csv(output_csv)
    assert sorted(df.columns.tolist()) == mmsplice_splicemap_cols



from conftest import outlier_dl_multi, outlier_model

def test_multi_sample_predict(outlier_dl_multi, outlier_model):

    results = outlier_model.predict_on_dataloader(outlier_dl_multi)
    print(results.df_mmsplice.shape)
    