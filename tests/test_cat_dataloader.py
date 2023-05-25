import pytest
from deepdiff import DeepDiff
from absplice import CatInference
from splicemap import SpliceCountTable as CountTable
from conftest import ref_table5_kn_testis, ref_table3_kn_testis, \
    ref_table5_kn_lung, ref_table3_kn_lung, \
    ref_table5_kn_blood, ref_table3_kn_blood, \
    count_cat_file_lymphocytes, count_cat_file_blood, \
    vcf_file, multi_vcf_file, fasta_file


@pytest.fixture
def cat_dl5():
    return CatInference(
        splicemap5=[ref_table5_kn_testis],
        count_cat=count_cat_file_lymphocytes,
    )


@pytest.fixture
def cat_dl3():
    return CatInference(
        splicemap3=[ref_table3_kn_testis],
        count_cat=count_cat_file_lymphocytes,
    )


# # only save necessary junctions
# # cat_dl[0].ct.df.loc[set.union(*cat_dl[0].common_junctions5).union(set.union(*cat_dl[0].common_junctions3))]
# def test_save_necessary_junctions():
    
#     count_cat_file_lymphocytes_complete = 'tests/data/create_test_data/full_data/backup/test_count_table_cat_chrom17_lymphocytes.csv'
#     count_cat_file_blood_complete = 'tests/data/create_test_data/full_data/backup/test_count_table_cat_chrom17_blood.csv'
    
#     cat_dl_complete = [
#         CatInference(splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
#                      splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
#                      count_cat=count_cat_file_lymphocytes_complete,
#                      name='lymphocytes'),
#         CatInference(splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
#                      splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
#                      count_cat=count_cat_file_blood_complete,
#                      name='blood'),
#     ]
    
#     cat_dl_complete[0].ct.df.loc[set.union(*cat_dl_complete[0].common_junctions5).union(set.union(*cat_dl_complete[0].common_junctions3))].to_csv(count_cat_file_lymphocytes, index=False)
#     cat_dl_complete[1].ct.df.loc[set.union(*cat_dl_complete[1].common_junctions5).union(set.union(*cat_dl_complete[1].common_junctions3))].to_csv(count_cat_file_blood, index=False)
    

def test_cat_dataloader_init(cat_dl):
    assert len(cat_dl) == 2
    assert [cat_dl[0].splicemaps5[i].method for i in range(len(cat_dl[0].splicemaps5))] == ['kn', 'kn']
    assert cat_dl[0].combined_splicemap5.shape[0] == 6
    assert sorted([cat_dl[0].splicemaps5[i].name for i in range(len(cat_dl[0].splicemaps5))]) \
        == sorted(['Testis', 'Lung'])

    assert [cat_dl[0].splicemaps3[i].method for i in range(len(cat_dl[0].splicemaps3))] == ['kn', 'kn']
    assert cat_dl[0].combined_splicemap3.shape[0] == 5
    assert sorted([cat_dl[0].splicemaps3[i].name for i in range(len(cat_dl[0].splicemaps3))]) \
        == sorted(['Testis', 'Lung'])


def test_cat_dataloader_common(cat_dl):
    assert len(set.union(*cat_dl[0].common_junctions5)) == 5
    assert len(set.union(*cat_dl[0].common_junctions3)) == 5


def test_cat_dataloader_common5(cat_dl5):
    assert len(cat_dl5.common_junctions5[0]) == 4


def test_cat_dataloader_common3(cat_dl3):
    assert len(cat_dl3.common_junctions3[0]) == 4


def test_cat_dataloader_sample_mapping():
    sample_mapping = {
        'NA00001': 'new_name1',
        'NA00002': 'new_name2',
        'NA00003': 'new_name3'
    }

    cat_dl_map = CatInference(
        splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
        splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
        count_cat=count_cat_file_lymphocytes,
        sample_mapping=sample_mapping)

    assert sorted(list(cat_dl_map.samples)) \
        == sorted(list({'new_name1', 'new_name2', 'new_name3'}))


def test_cat_dataloader_common(cat_dl):
    assert ('17:41201211-41203079:-', 'ENSG00000012048', 'Testis') \
        in set(cat_dl[1].common5.set_index(['junctions', 'gene_id', 'tissue']).index)
    
    assert not ('17:41201211-41203079:-', 'ENSG00000198496', 'Testis') \
        in set(cat_dl[1].common5.set_index(['junctions', 'gene_id', 'tissue']).index)
    
    assert not ('17:41201917-4120307900000:-', 'ENSG00000012048', 'Testis') \
        in set(cat_dl[0].common5.set_index(['junctions', 'gene_id', 'tissue']).index)
    
    # psi5
    for i in range(len(cat_dl)):
        df_common = cat_dl[i].common5
        for j in range(len(cat_dl[i].tissues5)):
            assert sorted(set(df_common[df_common['tissue'] == cat_dl[i].tissues5[j]]['junctions'])) \
                == sorted(set(cat_dl[i].common_junctions5[j]))
                
            _df = cat_dl[i].splicemaps5[j].df             
            assert sorted(set(df_common[df_common['tissue'] == cat_dl[i].tissues5[j]].set_index(['junctions', 'gene_id']).index)) \
                == sorted(set(_df[_df['junctions'].isin(cat_dl[i].common_junctions5[j])].set_index(['junctions', 'gene_id']).index))
       
    # psi3         
    for i in range(len(cat_dl)):
        df_common = cat_dl[i].common3
        for j in range(len(cat_dl[i].tissues3)):
            assert sorted(set(df_common[df_common['tissue'] == cat_dl[i].tissues3[j]]['junctions'])) \
                == sorted(set(cat_dl[i].common_junctions3[j]))
                
            _df = cat_dl[i].splicemaps3[j].df             
            assert sorted(set(df_common[df_common['tissue'] == cat_dl[i].tissues3[j]].set_index(['junctions', 'gene_id']).index)) \
                == sorted(set(_df[_df['junctions'].isin(cat_dl[i].common_junctions3[j])].set_index(['junctions', 'gene_id']).index))
    
                
def test_cat_dataloader_infer(cat_dl):
    # cat_dl[0].splicemaps5[0].df[['junctions', 'gene_id']]
    row = cat_dl[0].infer('17:41201211-41203079:-', 'ENSG00000012048', 'Testis', 'NA00002', 'psi5')
    assert row =={
        'junction': '17:41201211-41203079:-',
        'gene_id': 'ENSG00000012048',
        'sample': 'NA00002',
        'tissue': 'Testis',
        'tissue_cat': 'lymphocytes',
        'count_cat': 50,
        'delta_logit_psi_cat': -0.5108256237659907,
        'delta_psi_cat': -0.12500000000000006,
        'k_cat': 250,
        'median_n_cat': 300.0,
        'n_cat': 850,
        'psi_cat': 0.2,
        'ref_psi_cat': 0.29411764705882354}
    
    # ref_psi_cat is equal to ref_psi of target tissue -> delta_psi_cat is ref_psi_cat - psi_cat
    # cat_dl[0].splicemaps5[1].df[['junctions', 'gene_id']]
    row = cat_dl[0].infer('17:41277787-41290673:+', 'ENSG00000198496', 'Lung', 'NA00002', 'psi5')
    assert row =={
        'junction': '17:41277787-41290673:+',
        'gene_id': 'ENSG00000198496',
        'sample': 'NA00002',
        'tissue': 'Lung',
        'tissue_cat': 'lymphocytes',
        'count_cat': 20,
        'delta_logit_psi_cat': -1.6094379124341003,
        'delta_psi_cat': -0.3333333333333333,
        'k_cat': 300,
        'median_n_cat': 240.0,
        'n_cat': 600,
        'psi_cat': 0.16666666666666666,
        'ref_psi_cat': 0.5}
    
    # cat_dl[1].splicemaps5[0].df[['junctions', 'gene_id']]
    row = cat_dl[1].infer('17:41201211-41203079:-', 'ENSG00000012048', 'Testis', 'NA00002', 'psi5')
    assert row =={
        'junction': '17:41201211-41203079:-',
        'gene_id': 'ENSG00000012048',
        'sample': 'NA00002',
        'tissue': 'Testis',
        'tissue_cat': 'blood',
        'count_cat': 7,
        'delta_logit_psi_cat': 0.0,
        'delta_psi_cat': 0.0,
        'k_cat': 14,
        'median_n_cat': 5.0,
        'n_cat': 14,
        'psi_cat': 1.0,
        'ref_psi_cat': 1.0}


def test_cat_dataloader_contains(cat_dl):
    assert cat_dl[0].contains('NA00001')
    assert not cat_dl[0].contains('NA00005')
    

# def test_cat_dataloader_infer_only_one_cat(cat_dl):
#     junction_id = '17:41201917-41203079:-'
#     sample = 'NA00002'
#     tissue = 'Testis'
#     event_type = 'psi5'

#     row = cat_dl[0].infer(junction_id, sample, tissue, event_type)
#     assert row == {
#         'junction': '17:41201917-41203079:-',
#         'sample': 'NA00002',
#         'tissue': 'Testis',
#         'tissue_cat': 'lymphocytes',
#         'count_cat': 0,
#         'delta_logit_psi_cat': 0.0,
#         'delta_psi_cat': 1.734723475976807e-18,
#         'k_cat': 0,
#         'median_n_cat': 43.0,
#         'n_cat': 166,
#         'psi_cat': 0.0,
#         'ref_psi_cat': 0.0
#     }
#     assert cat_dl[1].contains(junction_id, sample, tissue, event_type) == False
    
def test_cat_dataloader_infer_splicemap_cat_with_splicemap_cat():
    cat_dl_splicemap_cat = CatInference(
        splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
        splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
        count_cat=count_cat_file_blood,
        splicemap_cat5=ref_table5_kn_blood,
        splicemap_cat3=ref_table3_kn_blood,
        name='blood',
    )
    
    cat_dl_no_splicemap_cat = CatInference(
        splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
        splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
        count_cat=count_cat_file_blood,
        name='blood',
    )
    
    # First junction is in count table and in splicemap (here median_n of higher stat power is used (from splicemap))
    junction_id = '17:41277787-41283224:+'
    gene_id = 'ENSG00000198496'
    # gene_id = 'ENSG00000012048' #not in
    sample = 'NA00002'
    tissue = 'Testis'
    event_type = 'psi5'

    assert junction_id in set(cat_dl_splicemap_cat.splicemap5_cat.df['junctions'])
    assert cat_dl_splicemap_cat.contains(sample) == True
    row = cat_dl_splicemap_cat.infer(junction_id, gene_id, tissue, sample, event_type)
    expected = {
        'junction': '17:41277787-41283224:+',
        'gene_id': 'ENSG00000198496',
        'sample': 'NA00002',
        'tissue': 'Testis',
        'tissue_cat': 'blood',
        'count_cat': 4,
        'delta_logit_psi_cat': -3.208825489014699,
        'delta_psi_cat': -0.18999999999999995,
        'k_cat': 12,
        'median_n_cat': 12.0,
        'n_cat': 14,
        'psi_cat': 0.8,
        'ref_psi_cat': 1.0
    }
    diff = DeepDiff(row, expected, significant_digits=10)
    assert not diff
    assert cat_dl_splicemap_cat.infer(junction_id, gene_id, tissue, sample, event_type) != \
        cat_dl_no_splicemap_cat.infer(junction_id, gene_id, tissue, sample, event_type)
    

def test_cat_dataloader_infer_splicemap_cat_with_splicemap_cat_not_in_splicemap():
    cat_dl_splicemap_cat = CatInference(
        splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
        splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
        count_cat=count_cat_file_blood,
        splicemap_cat5=ref_table5_kn_blood,
        splicemap_cat3=ref_table3_kn_blood,
        name='blood',
    )
    
    cat_dl_no_splicemap_cat = CatInference(
        splicemap5=[ref_table5_kn_testis, ref_table5_kn_lung],
        splicemap3=[ref_table3_kn_testis, ref_table3_kn_lung],
        count_cat=count_cat_file_blood,
        name='blood',
    )

    # junction_id not in SpliceMap of CAT, but in count table of CAT and in SpliceMap of target tissue -> use median_n from count table
    junction_id = '17:41201200-41205000:-'
    gene_id = 'ENSG00000012048'
    sample = 'NA00002'
    tissue = 'Testis'
    event_type = 'psi5'
    
    assert junction_id not in set(cat_dl_splicemap_cat.splicemap5_cat.df['junctions'])
    assert cat_dl_splicemap_cat.contains(sample) == True
    # _df = cat_dl_splicemap_cat.splicemaps5[cat_dl_splicemap_cat.tissues5.index(tissue)].df
    # _df[['junctions', 'gene_id']]
    row = cat_dl_splicemap_cat.infer(junction_id, gene_id, tissue, sample, event_type)
    assert row == {
        'junction': '17:41201200-41205000:-',
        'gene_id': 'ENSG00000012048',
        'sample': 'NA00002',
        'tissue': 'Testis',
        'tissue_cat': 'blood',
        'count_cat': 100,
        'delta_logit_psi_cat': 0.0,
        'delta_psi_cat': 0.0,
        'k_cat': 300,
        'median_n_cat': 100.0,
        'n_cat': 300,
        'psi_cat': 1.0,
        'ref_psi_cat': 1.0
    }
    assert cat_dl_splicemap_cat.infer(junction_id, gene_id, tissue, sample, event_type) == \
        cat_dl_no_splicemap_cat.infer(junction_id, gene_id, tissue, sample, event_type)


# def test_cat_dataloader_infer_all(cat_dl):

#     df = cat_dl[0].infer_all('psi5')

#     common_junctions = [i for sublist in cat_dl[0].common_junctions5 for i in sublist]
#     assert len(set(common_junctions)) == len(set(df.index.get_level_values('junction')))

#     tissues = cat_dl[0].tissues5
#     assert len(set(tissues)) == len(set(df.index.get_level_values('tissue')))

#     samples = cat_dl[0].samples
#     assert len(set(samples)) == len(set(df.index.get_level_values('sample')))

#     assert df.shape[0] <= len(set(common_junctions)) * len(set(tissues)) * len(set(samples))
