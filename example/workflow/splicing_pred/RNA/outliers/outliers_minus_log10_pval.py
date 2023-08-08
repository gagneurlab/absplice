import pandas as pd
import numpy as np
from absplice.utils import get_abs_max_rows

def minus_log10(x):
    return -np.log10(x)


df = pd.read_csv(snakemake.input['outlier_with_variant'])
df['tissue_cat'] = snakemake.wildcards['tissue_cat']

df['pValueGene_g_minus_log10'] = df['pValueGene_g'].apply(lambda x: minus_log10(x))

df = get_abs_max_rows(
    df.set_index(['variant', 'gene_id', 'sample']), 
    ['variant', 'gene_id', 'sample'], 
    'pValueGene_g_minus_log10'
).reset_index()

df = df[[
    'variant', 'gene_id', 'sample',
    'pValueGene_g_minus_log10', 'tissue_cat'
]]

df.to_csv(snakemake.output['outlier_cat_pval'], index=False)
    