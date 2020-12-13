import pandas as pd
import numpy as np


def get_abs_max_rows(df, groupby, max_col):
    df = df.reset_index()
    _df = df.copy()
    _df[max_col] = _df[max_col].abs()
    max_scores = _df.groupby(groupby)[max_col].idxmax()
    return df.iloc[max_scores.values].set_index(groupby)


def expit(x):
    return 1. / (1. + np.exp(-x))


def delta_logit_PSI_to_delta_PSI(delta_logit_psi, ref_psi,
                                 genotype=None, clip_threshold=0.001):
    ref_psi = clip(ref_psi, clip_threshold)
    pred_psi = expit(delta_logit_psi + logit(ref_psi))

    if genotype is not None:
        pred_psi = np.where(np.array(genotype) == 1,
                            (pred_psi + ref_psi) / 2,
                            pred_psi)

    return pred_psi - ref_psi


def clip(x, clip_threshold=0.01):
    return np.clip(x, clip_threshold, 1 - clip_threshold)


def logit(x, clip_threshold=0.01):
    x = clip(x, clip_threshold=clip_threshold)
    return np.log(x) - np.log(1 - x)
