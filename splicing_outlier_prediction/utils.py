import pandas as pd
import numpy as np
import pathlib


def get_abs_max_rows(df, groupby, max_col, dropna=True):
    df = df.reset_index()
    _df = df.copy()
    _df[max_col] = _df[max_col].abs()
    max_scores = _df.groupby(groupby, dropna=dropna)[max_col].idxmax()
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


def normalize_gene_annotation(df, gene_map, key='gene_name', value='gene_id'):
    if isinstance(gene_map, dict):
        pass
    elif isinstance(gene_map, pd.DataFrame):
        gene_map = dict(zip(gene_map[key], gene_map[value]))
    else:
        TypeError("gene_mapping needs to be dictionary of pandas DataFrame")
    df[value] = df[key].map(gene_map)
    return df


def read_csv(path, **kwargs):
    if not isinstance(path, pathlib.PosixPath):
        path = pathlib.Path(path)
    if path.suffix.lower() == '.csv':
        return pd.read_csv(path, **kwargs)
    elif path.suffix.lower() == '.tsv':
        return pd.read_csv(path, sep='\t', **kwargs)
    elif path.suffix.lower() == '.parquet':
        return pd.read_parquet(path, **kwargs)
    else:
        raise ValueError("unknown file ending.")
    
    
def filter_samples_with_RNA_seq(df, samples_for_tissue):
    """
        samples_for_tissue: Dict, keys: tissue, values: samples with RNA-seq for respective tissue
    """
    l = list()
    for tissue, samples in samples_for_tissue.items():
        df_with_RNA = df[(df['tissue'] == tissue) & (df['sample'].isin(samples))]
        l.append(df_with_RNA)
    return pd.concat(l)