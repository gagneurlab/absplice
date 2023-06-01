import pandas as pd
import pyranges as pr
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors import variants_to_pyranges
from absplice.utils import get_unique_splice_site_intervals_in_event, intervals_to_pyranges

# read in variants
df_vars = pd.read_csv(snakemake.input['variants'])

# read in outliers
df = pd.read_csv(snakemake.input['outliers_signif'])

df['junctions'] = df['junctions_j']
df['events'] = df['events_j']

df = df.set_index('junctions')
df['events'] = df['events'].str.split(';')

df['interval'] = df['events'].apply(lambda x: get_unique_splice_site_intervals_in_event(x, overhang=(250, 250)))
df_intervals_all = df.explode('interval')[['interval']]

# get unique intervals into pyranges
i_list = df_intervals_all['interval'].values
i_list = list(set(i_list))
pr_intervals_unique = intervals_to_pyranges(i_list)
df_intervals_unique = pr_intervals_unique.df

# get unique variants into pyranges
v_list = df_vars['variant'].apply(lambda x: Variant.from_str(x)).values
v_list = list(set(v_list))
pr_vars = variants_to_pyranges(v_list)

# get overlap of intervals and variants
pr_intervals_with_var = pr_intervals_unique.join(pr_vars)
df_intervals_with_var = pr_intervals_with_var.df

# interval to string for joining
df_intervals_all['interval'] = df_intervals_all['interval'].apply(lambda x: x.__str__())
df_intervals_unique['interval'] = df_intervals_unique['interval'].apply(lambda x: x.__str__())
df_intervals_with_var['interval'] = df_intervals_with_var['interval'].apply(lambda x: x.__str__())

# get junction information from df_intervals_all
df_intervals_all = df_intervals_all.reset_index().set_index('interval')
df_intervals_with_var = df_intervals_with_var.set_index('interval')
df_outlier_with_variant = df_intervals_all.join(
    df_intervals_with_var, how='inner').reset_index()[['junctions', 'variant']].drop_duplicates()

# get variant informations from df_vars
df_outlier_with_variant['variant'] = df_outlier_with_variant['variant'].apply(lambda x: x.__str__())
df_outlier_with_variant = df_outlier_with_variant.set_index('variant').join(
    df_vars.set_index('variant')).reset_index()

# add info from outliers
df = df.reset_index().drop(columns=['events_j', 'interval'])
df_outlier_with_variant = df_outlier_with_variant.set_index(['junctions', 'sample']).join(
    df.set_index(['junctions', 'sample']), how='inner').reset_index()

# save
df_outlier_with_variant.to_csv(snakemake.output['outlier_with_variant'], index=False)
