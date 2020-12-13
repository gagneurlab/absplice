import pandas as pd
from splicing_outlier_prediction.utils import get_abs_max_rows


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
