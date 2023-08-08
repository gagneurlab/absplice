import pandas as pd
from absplice.utils import annotate_junctions_DROP

# Read in sample mapping from RNA_ID to DNA_ID
df = pd.read_csv(snakemake.input['DROP_annotation'], sep='\t')
# subset for accessible tissue
df = df[
    df['DROP_GROUP'].apply(lambda x: x.split(','))\
        .apply(lambda x: snakemake.wildcards['tissue_cat'] in x)
]
sampleID_map_RNA_to_DNA = pd.Series(df['DNA_ID'].values,index=df['sampleID']).to_dict()

# gene level outliers
df_gene_level = pd.read_csv(snakemake.input['gene_level'], sep='\t')\
                    .rename(columns={'sampleID': 'sample'})
df_gene_level['sample'] = df_gene_level['sample'].astype(str)
df_gene_level['sample'] = df_gene_level['sample'].map(sampleID_map_RNA_to_DNA)
df_gene_level = annotate_junctions_DROP(df_gene_level)

# annotate gene_id
gene_map = pd.read_csv(snakemake.input['gene_map'], sep='\t')
gene_map = dict(zip(gene_map['gene_name'], gene_map['gene_id']))
df_gene_level['gene_id'] = df_gene_level['hgncSymbol'].map(gene_map)

# filter for protein coding genes
df_coding = pd.read_csv(snakemake.input['coding_genes'])
df_gene_level = df_gene_level[
    df_gene_level['gene_id'].isin(set(df_coding['gene_id']))
]

# filter by cutoffs
df_gene_level = df_gene_level[
    df_gene_level['padjustGene'] <= snakemake.params['padjustGene_cutoff']
]

# save results
df_gene_level.to_csv(snakemake.output['result'], index=False)
