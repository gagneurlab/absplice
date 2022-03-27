import pandas as pd
from splicing_outlier_prediction.utils import get_abs_max_rows, filter_samples_with_RNA_seq, read_spliceai_vcf, dtype_columns_spliceai
from splicing_outlier_prediction import SplicingOutlierResult
from conftest import gene_map, gene_tpm, spliceai_path, mmsplice_path, spliceai_vcf_path, spliceai_vcf_path2


def test_get_max_rows():
    df = pd.DataFrame({
        'junction': ['j1', 'j1', 'j2', 'j2', 'j3'],
        'sample': ['s1', 's1', 's2', 's2', 's2'],
        'score': [10, 20, -40, 30, 10]
    }).set_index('junction')

    pd.testing.assert_frame_equal(
        get_abs_max_rows(df, ['junction', 'sample'], 'score'),
        pd.DataFrame({
            'junction': ['j1', 'j2', 'j3'],
            'sample': ['s1', 's2', 's2'],
            'score': [20, -40, 10]
        }).set_index(['junction', 'sample'])
    )


# def test_outlier_results_filter_samples_with_RNA_seq(outlier_results_multi, outlier_model):
def test_outlier_results_filter_samples_with_RNA_seq(df_var_samples, outlier_model, gene_map, gene_tpm):
    samples_for_tissue = {
        'Testis': ['NA00002'],
        'Lung': ['NA00002', 'NA00003']
    }
    
    results = SplicingOutlierResult(
        df_mmsplice=mmsplice_path, 
        df_spliceai=spliceai_path, 
        gene_map=gene_map,
        gene_tpm=gene_tpm
    )

    results.add_samples(df_var_samples)

    assert results.df_mmsplice[['tissue', 'sample']].groupby('tissue')['sample'].apply(lambda x: ';'.join(sorted(list(set(x))))).to_dict() == \
        {
            # 'Lung': 'NA00002;NA00003', 
            'Testis': 'NA00002;NA00003'}
        
    df_spliceai_tpm = results._add_tissue_info_to_spliceai()
    assert df_spliceai_tpm[['tissue', 'sample']].groupby('tissue')['sample'].apply(lambda x: ';'.join(sorted(list(set(x))))).to_dict() == \
        {
            # 'Lung': 'NA00002;NA00003', 
            'Testis': 'NA00002;NA00003'}

    df_mmsplice = filter_samples_with_RNA_seq(results.df_mmsplice, samples_for_tissue)
    df_spliceai = filter_samples_with_RNA_seq(df_spliceai_tpm, samples_for_tissue)

    assert df_mmsplice[['tissue', 'sample']].groupby('tissue')['sample'].apply(lambda x: ';'.join(sorted(list(set(x))))).to_dict() == \
        {
            # 'Lung': 'NA00002;NA00003', 
            'Testis': 'NA00002'}
    assert df_spliceai[['tissue', 'sample']].groupby('tissue')['sample'].apply(lambda x: ';'.join(sorted(list(set(x))))).to_dict() == \
        {
            # 'Lung': 'NA00002;NA00003', 
            'Testis': 'NA00002'}
        
    
def test_utils_read_spliceai_vcf():
    df = read_spliceai_vcf(spliceai_vcf_path2)
    df_compare = pd.DataFrame({
        'variant': ['17:41201201:TTC>CA', '17:41276032:T>A', '17:41279042:A>GA'],
        'gene_name': ['OR4F5', 'EXOSC3', 'test'],
        'delta_score': [0.01, 0.00, 0.00],
        'acceptor_gain': [0.01, 0.00, 0.00],
        'acceptor_loss': [0.00, 0.00, 0.00],
        'donor_gain': [0.00, 0.00, 0.00],
        'donor_loss': [0.00, 0.00, 0.00],
        'acceptor_gain_position': [42, 0, 0],
        'acceptor_loss_positiin': [25, -13, -13],
        'donor_gain_position': [24, -44, -44],
        'donor_loss_position': [2, -12, -12]
    })
    for col in df_compare.columns:
        if col in dtype_columns_spliceai.keys():
            df_compare = df_compare.astype({col: dtype_columns_spliceai[col]})
    pd.testing.assert_frame_equal(df, df_compare)