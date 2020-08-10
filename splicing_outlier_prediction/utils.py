import pandas as pd


def get_abs_max_rows(df, groupby, max_col):
    df = df.copy()
    df[max_col] = df[max_col].abs()
    df = df.reset_index()
    max_scores = df.groupby(groupby)[max_col].idxmax()
    return df.iloc[max_scores.values].set_index(groupby)
