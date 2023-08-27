import pandas as pd
import numpy as np
from pyranges import PyRanges
import pathlib
from kipoiseq import Interval
from kipoiseq.extractors.vcf import MultiSampleVCF
from kipoiseq.extractors.vcf_query import BaseVariantQuery
from collections import namedtuple
from typing import List


# def get_abs_max_rows(df, groupby, max_col, dropna=True):
#     if len(max_col) == 1:
#         return df.reset_index().sort_values(by=max_col, key=abs, ascending=False).drop_duplicates(subset=groupby).set_index(groupby)
#     else:
#         return df.reset_index() \
#             .sort_values(by=max_col, key=abs, ascending=[False, False]) \
#             .drop_duplicates(subset=groupby) \
#             .set_index(groupby)

def get_abs_max_rows(df, groupby, max_col, dropna=True):

    if type(max_col) == str:
        return df.reset_index() \
            .sort_values(by=max_col, key=abs, ascending=False) \
            .drop_duplicates(subset=groupby) \
            .set_index(groupby)
    else:
        return df.reset_index() \
            .sort_values(by=max_col, key=abs, ascending=[False, False]) \
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
    return pd.concat([df, new_row])


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
                results = [0 if e == '.' else e for e in results]
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


def _add_variant_col(row):
    if '#Chrom' in row:
        return str(row['#Chrom']).replace('chr', '') + ':' + str(row['Pos']) + ':' + row['Ref'] + '>' + row['Alt']
    elif 'chrom' in row:
        return str(row['chrom']).replace('chr', '') + ':' + str(row['pos']) + ':' + row['ref'] + '>' + row['alt']
    else:
        raise NotImplementedError()


def _add_variant(df):
    if 'variant' in df.reset_index().columns:
        return df
    else:
        df['variant'] = df.apply(lambda x: _add_variant_col(x), axis=1)
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


def read_absplice(path, **kwargs):
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
    return df


def annotate_junctions_DROP(df):
    if 'chr' not in df['seqnames'].astype(str).values[0]:
        df['seqnames'] = df['seqnames'].astype(str)
        df['seqnames'] = df['seqnames'].apply(lambda x: 'chr' + x)
    df['junctions'] = df['seqnames'].astype(str) + ':' \
            + (df['start'] - 1).astype(str) + '-' \
            + df['end'].astype(str) + ':' \
            + df['strand'].map(lambda x: x if x == '-' else '+')
    return df


class VariantMafFilter(BaseVariantQuery):

    def __init__(self, cutoff, population=None, not_in_population=True):
        self.population = population
        self.cutoff = cutoff
        self.not_in_population = not_in_population

    def __call__(self, v):
        if self.population:
            v = str(v)
            return self.population[v] <= self.cutoff \
                if v in self.population else self.not_in_population
        else:
            return v.source.aaf <= self.cutoff


class PrivateVariantFilter(BaseVariantQuery):

    def __init__(self, vcf, max_num_samples=2):
        self.vcf = vcf
        self.max_num_sample = max_num_samples

    def __call__(self, v):
        return len(self.vcf.get_samples(v)) <= self.max_num_sample


class ReadDepthFilter(BaseVariantQuery):

    def __init__(self, vcf, min_read=10, sample_id=None):
        self.sample_mapping = vcf.sample_mapping
        self.min_read = min_read
        self.sample_id = sample_id

    def __call__(self, v):
        gt_depth = v.source.gt_alt_depths
        if self.sample_id:
            depth = gt_depth[self.sample_mapping[self.sample_id]]
        else:
            depth = min(gt_depth)
        return depth >= self.min_read


class GQFilter(BaseVariantQuery):

    def __init__(self, vcf, min_GQ=80, sample_id=None):
        self.sample_mapping = vcf.sample_mapping
        self.min_GQ = min_GQ
        self.sample_id = sample_id

    def __call__(self, v):
        gt_quals = v.source.gt_quals
        if self.sample_id:
            GQ = gt_quals[self.sample_mapping[self.sample_id]]
        else:
            GQ = max(gt_quals)
        return GQ >= self.min_GQ
    
    
class LongVariantFilter(BaseVariantQuery):

    def __init__(self, vcf, max_length=10):
        self.vcf = vcf
        self.max_length = max_length

    def __call__(self, v):
        v = left_normalized(v)
        length = max(len(v.ref), len(v.alt))
        return length < self.max_length
    
    
class Junction(Interval):

    @property
    def acceptor(self):
        return self.start if self.strand == '-' else self.end

    @property
    def donor(self):
        return self.end if self.strand == '-' else self.start

    def dinucleotide_region(self):
        return Interval(self.chrom, self.start, self.start + 2), \
            Interval(self.chrom, self.end - 2, self.end)

    def acceptor_region(self, overhang=(250, 250)):
        return Interval(self.chrom, self.acceptor,
                        self.acceptor, strand=self.strand) \
            .slop(upstream=overhang[0], downstream=overhang[1])

    def donor_region(self, overhang=(250, 250)):
        return Interval(self.chrom, self.donor,
                        self.donor, strand=self.strand) \
            .slop(upstream=overhang[0], downstream=overhang[1])
            
    
def get_splice_site_intervals(junction, overhang=(250, 250)):
    junction = Junction.from_str(junction) if type(
        junction) == str else junction

    acceptor = junction.acceptor_region(overhang=overhang)
    donor = junction.donor_region(overhang=overhang)
    return [acceptor, donor]

def get_unique_splice_site_intervals_in_event(event, overhang=(250, 250)):
    sites = list()
    for junction in event:
        sites.append(get_splice_site_intervals(junction, overhang))
    sites = [item for sublist in sites for item in sublist]
    sites = list(set(sites))
    return sites
    
def intervals_to_pyranges(intervals: List[Interval]) -> PyRanges:
    """
    Create pyrange object given list of intervals objects.
    Args:
      intervals: list of interval objects have CHROM, START, END, properties.
    """
    import pyranges
    df = pd.DataFrame([
        (
            i.chrom,
            i.start,
            i.end,
            i
        )
        for i in intervals
    ], columns=['Chromosome', 'Start', 'End', 'interval'])
    return pyranges.PyRanges(df)
