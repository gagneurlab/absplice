import pytest
from splicing_outlier_prediction import CatSpliceOutlier


@pytest.fixture
def cat_outlier_model():
    return CatSpliceOutlier()


def test_cat_splicing_outlier_on_batch(cat_outlier_model, cat_dl):
    batch = next(cat_dl.batch_iter())
    df = cat_outlier_model.predict_on_batch(batch, cat_dl)

    assert df.columns.tolist() == [
        'junction', 'Chromosome', 'Start', 'End', 'Strand', 'events',
        'splice_site', 'ref_psi', 'k', 'n', 'gene_id', 'gene_name', 'weak',
        'transcript_id', 'gene_type', 'cat_ref_psi', 'cat_k', 'cat_n',
        'NA00001_cat_psi', 'NA00001_cat_delta_logit_psi',
        'NA00001_target_delta_psi', 'NA00002_cat_psi',
        'NA00002_cat_delta_logit_psi', 'NA00002_target_delta_psi',
        'NA00003_cat_psi', 'NA00003_cat_delta_logit_psi',
        'NA00003_target_delta_psi'
    ]
