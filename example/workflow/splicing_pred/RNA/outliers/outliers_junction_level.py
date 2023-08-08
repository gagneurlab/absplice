import pandas as pd
from splicemap.count_table import SpliceCountTable as CountTable
from absplice.utils import annotate_junctions_DROP

# Read in sample mapping from RNA_ID to DNA_ID
df = pd.read_csv(snakemake.input['DROP_annotation'], sep='\t')
# subset for accessible tissue
df = df[
    df['DROP_GROUP'].apply(lambda x: x.split(','))\
        .apply(lambda x: snakemake.wildcards['tissue_cat'] in x)
]
sampleID_map_RNA_to_DNA = pd.Series(df['DNA_ID'].values,index=df['sampleID']).to_dict()

# junction level outliers
df_junction_level = pd.read_csv(snakemake.input['junction_level'], sep='\t')\
                    .rename(columns={'sampleID': 'sample'})
df_junction_level['sample'] = df_junction_level['sample'].astype(str)
df_junction_level['sample'] = df_junction_level['sample'].map(sampleID_map_RNA_to_DNA)
df_junction_level = annotate_junctions_DROP(df_junction_level)

# annotate gene_id
gene_map = pd.read_csv(snakemake.input['gene_map'], sep='\t')
gene_map = dict(zip(gene_map['gene_name'], gene_map['gene_id']))
df_junction_level['gene_id'] = df_junction_level['hgncSymbol'].map(gene_map)

# filter for protein coding genes
df_coding = pd.read_csv(snakemake.input['coding_genes'])
df_junction_level = df_junction_level[
    df_junction_level['gene_id'].isin(set(df_coding['gene_id']))
]

# filter by cutoffs
df_junction_level = df_junction_level[
    (df_junction_level['padjust'] <= snakemake.params['padjustJunction_cutoff'])
    & (df_junction_level['totalCounts'] >= snakemake.params['totalCounts_cutoff'])
    & (df_junction_level['deltaPsi'].abs() >= snakemake.params['delta_psi_cutoff'])
]

# add splice site and event information from count table
ct = CountTable.read_csv(snakemake.input['count_table_updated'])

df_junction_psi5 = df_junction_level[df_junction_level['type'] == 'psi5'].set_index('junctions')
df_junction_psi3 = df_junction_level[df_junction_level['type'] == 'psi3'].set_index('junctions')

df_junction_psi5 = df_junction_psi5.join(ct.splice_site5).join(ct.event5).reset_index()
df_junction_psi3 = df_junction_psi3.join(ct.splice_site3).join(ct.event3).reset_index()

df_junction_theta = df_junction_level[df_junction_level['type'] == 'theta']
df_junction_theta['events'] = df_junction_theta['junctions']

pd.concat([
    df_junction_psi5, 
    df_junction_psi3,
    df_junction_theta
]).to_csv(snakemake.output['result'], index=False)
