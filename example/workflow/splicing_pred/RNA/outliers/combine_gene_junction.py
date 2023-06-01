import pandas as pd

df_junction = pd.read_csv(snakemake.input['junction'])
df_gene = pd.read_csv(snakemake.input['gene'])

df_g_j = df_junction.set_index(['gene_id', 'sample']).add_suffix('_j')\
        .join(df_gene.set_index(['gene_id', 'sample']).add_suffix('_g'), how='inner')\
                .reset_index()

df_g_j.to_csv(snakemake.output['gene_junction'], index=False)
