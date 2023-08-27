import pandas as pd
from absplice.utils import get_abs_max_rows, filter_samples_with_RNA_seq, \
    read_cadd_splice, read_absplice, read_spliceai_vcf, dtype_columns_spliceai
from absplice import SplicingOutlierResult
from conftest import gene_map, gene_tpm, spliceai_path, mmsplice_path, \
    spliceai_vcf_path, spliceai_vcf_path2, cadd_splice_path#, \
        # absplice_precomputed_path


def test_get_max_rows():
    df_1 = pd.DataFrame({
        'junction': ['j1', 'j1', 'j2', 'j2', 'j3'],
        'sample': ['s1', 's1', 's2', 's2', 's2'],
        'score': [10, 20, -40, 30, 10]
    }).set_index('junction')

    pd.testing.assert_frame_equal(
        get_abs_max_rows(df_1, ['junction', 'sample'], 'score'),
        pd.DataFrame({
            'junction': ['j1', 'j2', 'j3'],
            'sample': ['s1', 's2', 's2'],
            'score': [20, -40, 10]
        }).sort_values(by='score', key=abs, ascending=False).set_index(['junction', 'sample'])
    )

    df_2 = pd.DataFrame({
        'junction': ['j1', 'j1', 'j2', 'j2', 'j3'],
        'sample': ['s1', 's1', 's2', 's2', 's2'],
        'score': [10, 20, -40, 30, -40], 
        'count': [15, 40, 5, 100, 150]
    }).set_index('junction')

    pd.testing.assert_frame_equal(
        get_abs_max_rows(df_2, ['junction', 'sample'], ['score', 'count']),
        pd.DataFrame({
            'junction': ['j3', 'j2', 'j1'],
            'sample': ['s2', 's2', 's1'],
            'score': [-40, -40, 20],
            'count': [150, 5, 40]
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
        'variant': ['17:41201201:TTC>CA', '17:41201201:TTC>CA', '17:41276032:T>A', '17:41279042:A>GA'],
        'gene_name': ['OR4F5', 'gene2', 'EXOSC3', 'test'],
        'delta_score': [0.01, 0.01, 0.00, 0.00],
        'acceptor_gain': [0.01, 0.01, 0.00, 0.00],
        'acceptor_loss': [0.00, 0.00, 0.00, 0.00],
        'donor_gain': [0.00, 0.00, 0.00, 0.00],
        'donor_loss': [0.00, 0.00, 0.00, 0.00],
        'acceptor_gain_position': [42, 42, 0, 0],
        'acceptor_loss_positiin': [25, 25, -13, -13],
        'donor_gain_position': [24, 24, -44, -44],
        'donor_loss_position': [2, 2, -12, -12]
    })
    for col in df_compare.columns:
        if col in dtype_columns_spliceai.keys():
            df_compare = df_compare.astype({col: dtype_columns_spliceai[col]})
    pd.testing.assert_frame_equal(df, df_compare)
    
    
def test_utils_read_cadd_splice():
    df = read_cadd_splice(cadd_splice_path, skiprows=1)
    df_compare = pd.DataFrame({
        '#Chrom': [17, 17, 17, 17, 17],
        'Pos': [41229271, 41236605, 41249878, 41271389, 41238435],
        'Ref': ['T', 'G', 'G', 'A', 'C'],
        'Alt': ['A', 'A', 'A', 'G', 'G'],
        'RawScore': [ 0.087061, -0.546492, -0.987847,  0.071266,  0.18075 ],
        'PHRED': [1.978, 0.095, 0.008, 1.846, 2.907],
        'gene_id': ['ENSG00000204873', 'ENSG00000204873', 'ENSG00000241595','ENSG00000212659', 'ENSG00000204873'],
        'variant': ['17:41229271:T>A', '17:41236605:G>A', '17:41249878:G>A', '17:41271389:A>G', '17:41238435:C>G'],
    })
    pd.testing.assert_frame_equal(df, df_compare)
    
    
# def test_utils_read_absplice():
#     df = read_absplice(absplice_precomputed_path)
#     df_compare = pd.DataFrame({
#         'gene_id': 'ENSG00000012048',
#         'tissue': 'Testis',
#         'delta_logit_psi': 0.0006480262366686,
#         'delta_psi': 6.413423046081057e-06,
#         'delta_score': 0.8,
#         'splice_site_is_expressed': 1,
#         'median_n': 51,
#         'AbSplice_DNA': 0.0305930295000895,
#         'chrom': 'chr17',
#         'pos': 41276032,
#         'ref': 'T',
#         'alt': 'A',
#         'start': 41276032,
#         'end': 41276033,
#         'variant': '17:41276032:T>A'
#     }, index=[0])
#     pd.testing.assert_frame_equal(df, df_compare)
    
