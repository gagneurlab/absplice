import pytest
import pandas as pd
from splicing_outlier_prediction.spliceAI import SpliceAI
from conftest import fasta_file, spliceai_db_path, multi_vcf_file


@pytest.fixture
def spliceai():
    return SpliceAI(fasta_file, 'grch37')


def test_SpliceAI_parse():
    score = SpliceAI.parse('CA|TTN|0.07|1.00|0.00|0.00|-7|-1|35|-29')

    assert score == SpliceAI.Score(
        'TTN', 1., 0.07, 1., 0., 0., -7, -1, 35, -29)


def test_SpliceAI_predict(spliceai):
    scores = spliceai.predict('17:34149615:A>T')
    # 'TAF15|0.25|0|0|0|11|29|11|33'

    assert scores[0] == spliceai.Score(
        gene_name='TAF15', delta_score=0.25,
        acceptor_gain=0.25, acceptor_loss=0.0,
        donor_gain=0.0, donor_loss=0.0,
        acceptor_gain_position=11.0, acceptor_loss_positiin=29.0,
        donor_gain_position=11.0, donor_loss_position=33.0)

    scores = spliceai.predict('17:883700:A>T')
    assert scores == []


def test_SpliceAI_predict_df(spliceai):
    df = spliceai.predict_df([
        '17:34149615:A>T',
        '17:76202989:CAGATTATCTTGAA>C',
        '17:883700:A>T'
    ])
    pd.testing.assert_frame_equal(df, pd.DataFrame({
        'variant': ['17:34149615:A>T', '17:76202989:CAGATTATCTTGAA>C'],
        'gene_name': ['TAF15', 'AFMID'],
        'delta_score': [0.25, 0.98],
        'acceptor_gain': [0.25, 0.74],
        'acceptor_loss': [0.0, 0.98],
        'donor_gain': [0.0, 0.0],
        'donor_loss': [0.0, 0.0],
        'acceptor_gain_position': [11.0, 27.0],
        'acceptor_loss_positiin': [29.0, 3.0],
        'donor_gain_position': [11.0, 30.0],
        'donor_loss_position': [33.0, 3.0]
    }).set_index('variant'))


@pytest.fixture
def spliceai_db():
    return SpliceAI(fasta_file, 'grch37', db_path=spliceai_db_path)


def test_SpliceAI_predict_with_db(spliceai_db, mocker):
    scores = spliceai_db.predict('17:34149615:A>T')
    assert scores[0] == spliceai_db.Score(
        gene_name='TAF15', delta_score=0.25,
        acceptor_gain=0.25, acceptor_loss=0.25,
        donor_gain=0.0, donor_loss=0.0,
        acceptor_gain_position=11.0, acceptor_loss_positiin=29.0,
        donor_gain_position=11.0, donor_loss_position=33.0)


def test_SpliceAI_predict_on_vcf(spliceai_db, tmp_path):
    spliceai_db.samples = True
    spliceai_db.quality = True
    output_csv = tmp_path / 'output.csv'
    spliceai_db.predict_save(multi_vcf_file, output_csv)
    df = pd.read_csv(output_csv)

    assert df.columns.tolist() == [
        'variant', 'gene_name', 'delta_score', 'acceptor_gain',
        'acceptor_loss', 'donor_gain', 'donor_loss', 'acceptor_gain_position',
        'acceptor_loss_positiin', 'donor_gain_position', 'donor_loss_position',
        'samples', 'GQ', 'DP_ALT'
    ]


def test_SpliceAI_predict_db_only():
    spliceai = SpliceAI(db_path=spliceai_db_path)
    df = spliceai.predict_df(['17:34149615:A>T', '17:34149617:A>T'])
    assert df.shape == (1, 10)

    spliceai = SpliceAI(db_path=spliceai_db_path)
    df = spliceai.predict_df(['chr17:34149615:A>T', 'chr17:34149617:A>T'])
    assert df.shape == (1, 10)
