import pandas as pd


def get_abs_max_rows(df, groupby, max_col):
    df = df.copy()
    df[max_col] = df[max_col].abs()
    df = df.reset_index()
    max_scores = df.groupby(groupby)[max_col].idxmax()
    return df.iloc[max_scores.values].set_index(groupby)


def delta_logit_PSI_to_delta_PSI(delta_logit_psi, ref_psi,
                                 genotype=None, clip_threshold=0.001):
    ref_psi = clip(ref_psi, clip_threshold)
    pred_psi = expit(delta_logit_psi + logit(ref_psi))

    if genotype is not None:
        pred_psi = np.where(np.array(genotype) == 1,
                            (pred_psi + ref_psi) / 2,
                            pred_psi)

    return pred_psi - ref_psi


def logit(x, clip_threshold=0.00001):
    x = clip(x, clip_threshold=clip_threshold)
    return np.log(x) - np.log(1 - x)