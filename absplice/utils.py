import pandas as pd
import numpy as np
import pathlib
from kipoiseq.extractors.vcf import MultiSampleVCF
from collections import namedtuple


def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)
    
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
    elif isinstance(gene_map, pathlib.PosixPath) or isinstance(gene_map, str):
        gene_map = read_csv(gene_map)
        gene_map = dict(zip(gene_map[key], gene_map[value]))
    elif isinstance(gene_map, pd.DataFrame):
        gene_map = dict(zip(gene_map[key], gene_map[value]))
    else:
        TypeError("gene_mapping needs to be dictionary, pandas DataFrame or path")
    df[value] = df[key].map(gene_map)
    return df


def read_csv(path, **kwargs):
    if isinstance(path, pd.DataFrame):
        return path
    else:
        if not isinstance(path, pathlib.PosixPath):
            path = pathlib.Path(path)
        if path.suffix.lower() == '.csv' or str(path).endswith('.csv.gz'):
            return pd.read_csv(path, **kwargs)
        elif path.suffix.lower() == '.tsv' or str(path).endswith('.tsv.gz'):
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


def inject_new_row(df, new_row_dict):
    new_row = pd.DataFrame(df[-1:].values, columns=df.columns)
    for k,v in new_row_dict.items():
        new_row[k] = v
    return df.append(new_row)


def read_spliceai(path, **kwargs):
    if isinstance(path, pd.DataFrame):
        return path
    else:
        if not isinstance(path, pathlib.PosixPath):
            path = pathlib.Path(path)
        if path.suffix.lower() == '.csv' or str(path).endswith('.csv.gz'):
            return pd.read_csv(path, **kwargs)
        elif path.suffix.lower() == '.tsv' or str(path).endswith('.tsv.gz'):
            return pd.read_csv(path, sep='\t', **kwargs)
        elif path.suffix.lower() == '.parquet':
            return pd.read_parquet(path, **kwargs)
        elif path.suffix.lower() == '.vcf' or str(path).endswith('.vcf.gz'):
            return read_spliceai_vcf(path)
        else:
            raise ValueError("unknown file ending.")


dtype_columns_spliceai = {
    'variant': pd.StringDtype(), 
    'gene_name': pd.StringDtype(), 
    'delta_score': 'float64',
    'acceptor_gain': 'float64', 
    'acceptor_loss': 'float64',
    'donor_gain': 'float64', 
    'donor_loss': 'float64',
    'acceptor_gain_position': 'Int64',
    'acceptor_loss_positiin': 'Int64',
    'donor_gain_position': 'Int64',
    'donor_loss_position': 'Int64',
}

def read_spliceai_vcf(path):
    columns = ['gene_name', 'delta_score',
                'acceptor_gain', 'acceptor_loss',
                'donor_gain', 'donor_loss',
                'acceptor_gain_position',
                'acceptor_loss_positiin',
                'donor_gain_position',
                'donor_loss_position']
    rows = list()
    for variant in MultiSampleVCF(path):
        row_all = variant.source.INFO.get('SpliceAI')
        if row_all:
            for row in row_all.split(','):
                results = row.split('|')[1:]
                scores = np.array(list(map(float, results[1:])))
                spliceai_info = [results[0], scores[:4].max(), *scores]
                rows.append({
                    **{'variant': str(variant)}, 
                    **dict(zip(columns, spliceai_info))})
    df = pd.DataFrame(rows)
        
    for col in df.columns:
        if col in dtype_columns_spliceai.keys():
            df = df.astype({col: dtype_columns_spliceai[col]})
            
    return df


def _add_variant_row(row):
    return str(row['#Chrom']) + ':' + str(row['Pos']) + ':' + row['Ref'] + '>' + row['Alt']


def _add_variant(df):
    if 'variant' in df.columns:
        return df
    else:
        df['variant'] = df.apply(lambda x: _add_variant_row(x), axis=1)
        # df['variant'] = df['variant'].astype(pd.StringDtype())
        return df
    
    
def _check_gene_id(df):
    if 'GeneID' in df.columns:
        df = df.rename(columns={'GeneID': 'gene_id'})
    if 'gene_id' not in df.columns:
        raise ValueError("Please add gene_id.")
    return df

        
def read_cadd_splice(path, **kwargs):
    if isinstance(path, pd.DataFrame):
        df = path
    else:
        if not isinstance(path, pathlib.PosixPath):
            path = pathlib.Path(path)
        if path.suffix.lower() == '.csv' or str(path).endswith('.csv.gz'):
            df = pd.read_csv(path, **kwargs)
        elif path.suffix.lower() == '.tsv' or str(path).endswith('.tsv.gz'):
            df = pd.read_csv(path, sep='\t', **kwargs)
        elif path.suffix.lower() == '.parquet':
            df  = pd.read_parquet(path, **kwargs)
        else:
            raise ValueError("unknown file ending.")
    df = _add_variant(df)
    df = _check_gene_id(df)
    return df
