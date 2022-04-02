# import random
# import pytest
# import pandas as pd
# import numpy as np
# from absplice import SpliceOutlier, SpliceOutlierDataloader, CatInference
# from absplice.ensemble import train_model_ebm
# from conftest import fasta_file, multi_vcf_file, \
#     ref_table5_kn_testis, ref_table3_kn_testis,  \
#     ref_table5_kn_lung, ref_table3_kn_lung, \
#     combined_ref_tables5_testis_lung, combined_ref_tables3_testis_lung, \
#     count_cat_file_lymphocytes,  count_cat_file_blood, \
#     spliceAI, pickle_DNA, pickle_DNA_CAT


# def test_splicing_outlier_result_predict_ensemble_DNA(outlier_results_multi):

#     results = outlier_results_multi
#     results.add_spliceAI(spliceAI)

#     features_DNA = ['delta_psi', 'delta_score', 'median_n', 'ref_psi']
#     results.predict_ensemble(pickle_DNA, results.gene, features_DNA)
#     assert sorted(results._ensemble.columns.tolist()) == sorted([
#         'junction', 'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
#         'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#         'gene_id', 'transcript_id', 'gene_type',  
#         'delta_psi', 'delta_logit_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#         'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#         'index', 'variant_spliceAI', 'delta_score', 'acceptor_gain',
#         'acceptor_loss', 'donor_gain', 'donor_loss', 'acceptor_gain_position',
#         'acceptor_loss_positiin', 'donor_gain_position', 'donor_loss_position',
#         'GQ', 'DP_ALT', 'ensemble_pred'
#     ])

#     assert results._ensemble['ensemble_pred'] is not None


# def test_splicing_outlier_result_predict_ensemble_DNA_CAT(outlier_results_multi, cat_dl):
    
#     results = outlier_results_multi
#     results.add_spliceAI(spliceAI)
#     results.infer_cat(cat_dl)

#     features_DNA_CAT = ['delta_psi', 'delta_psi_cat', 'delta_score',
#                         'ref_psi', 'median_n', 'psi_cat', 'ref_psi_cat']
#     results.predict_ensemble(
#         pickle_DNA_CAT, results.gene_cat_concat, features_DNA_CAT)
#     assert sorted(results._ensemble.columns.tolist()) == sorted([
#         'junction', 'event_type', 'variant', 'Chromosome', 'Start', 'End', 'Strand',
#         'events', 'splice_site', 'ref_psi', 'k', 'n', 'median_n',
#         'novel_junction', 'weak_site_acceptor', 'weak_site_donor',
#         'gene_id', 'transcript_id', 'gene_type',  
#         'delta_psi', 'delta_logit_psi',
#         'ref_acceptorIntron', 'ref_acceptor', 'ref_exon', 'ref_donor', 'ref_donorIntron',
#         'alt_acceptorIntron', 'alt_acceptor', 'alt_exon', 'alt_donor', 'alt_donorIntron',
#         'count_cat', 'delta_logit_psi_cat', 'delta_psi_cat', 'k_cat', 'median_n_cat', 'n_cat', 'psi_cat', 'ref_psi_cat', 'tissue_cat',
#         'index', 'variant_spliceAI', 'delta_score', 'acceptor_gain',
#         'acceptor_loss', 'donor_gain', 'donor_loss', 'acceptor_gain_position',
#         'acceptor_loss_positiin', 'donor_gain_position', 'donor_loss_position',
#         'GQ', 'DP_ALT', 'ensemble_pred'
#     ])

#     assert results._ensemble['ensemble_pred'] is not None


# def test_splicing_outlier_result_train_ensemble_DNA(outlier_results_multi):

#     results = outlier_results_multi
#     results.add_spliceAI(spliceAI)

#     features_DNA = ['delta_psi', 'delta_score', 'median_n', 'ref_psi']

#     results.gene
#     # results._gene = results._gene.fillna(0)
#     results.gene['outlier'] = np.random.randint(0, 2, results.gene.shape[0])

#     results_ensemble, models = train_model_ebm(
#         results.gene, features_DNA, feature_to_filter_na=None, nsplits=2)

#     assert sorted(results_ensemble.columns.tolist()) == sorted([
#         *results._gene.index.names, *features_DNA,
#         'fold', 'y_test', 'y_pred'
#     ])


# def test_splicing_outlier_result_train_ensemble_DNA_CAT(outlier_results_multi, cat_dl):
    
#     results = outlier_results_multi
#     results.add_spliceAI(spliceAI)
#     results.infer_cat(cat_dl)

#     features_DNA_CAT = ['delta_psi', 'delta_psi_cat', 'delta_score',
#                         'ref_psi', 'median_n', 'psi_cat', 'ref_psi_cat']

#     # results._gene = results._gene.fillna(0)
#     results.gene_cat_concat
#     results.gene_cat_concat['outlier'] = np.random.randint(
#         0, 2, results.gene_cat_concat.shape[0])

#     results_ensemble, models = train_model_ebm(
#         results.gene_cat_concat, features_DNA_CAT, feature_to_filter_na=None, nsplits=2)

#     assert sorted(results_ensemble.columns.tolist()) == sorted([
#         *results.gene.index.names, *features_DNA_CAT,
#         'fold', 'y_test', 'y_pred'
#     ])


# def test_splicing_outlier_result_train_ensemble_DNA_CAT_cross_apply(outlier_results_multi, cat_dl):
    
#     results = outlier_results_multi
#     results.add_spliceAI(spliceAI)
#     results.infer_cat(cat_dl)

#     features_DNA = ['delta_psi', 'delta_score', 'ref_psi', 'median_n']
#     features_blood = ['delta_psi_blood', 'k_blood',
#                       'median_n_blood', 'n_blood', 'psi_blood', 'ref_psi_blood']
#     features_lympho = ['delta_psi_lymphocytes', 'k_lymphocytes',
#                        'median_n_lymphocytes', 'n_lymphocytes', 'psi_lymphocytes', 'ref_psi_lymphocytes']
#     features_DNA_CAT = [*features_DNA, *features_blood, *features_lympho]
#     features_DNA_CAT_train = [*features_DNA, *features_blood]
#     features_DNA_CAT_test = [*features_DNA, *features_lympho]

#     # results._gene = results._gene.fillna(0)
#     results.gene_cat_features
#     results.gene_cat_features['outlier'] = np.random.randint(
#         0, 2, results.gene_cat_features.shape[0])

#     results_ensemble, models = train_model_ebm(results.gene_cat_features, features_DNA_CAT,
#                                                feature_to_filter_na=None, nsplits=2,
#                                                features_train=features_DNA_CAT_train, features_test=features_DNA_CAT_test)

#     assert sorted(results_ensemble.columns.tolist()) == sorted([
#         *results.gene.index.names, *features_DNA_CAT,
#         'fold', 'y_test', 'y_pred', 'y_pred_on_train_features'
#     ])
