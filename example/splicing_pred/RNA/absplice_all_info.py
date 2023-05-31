import pandas as pd

df_absplice_dna = pd.read_csv(snakemake.input['absplice_dna'])
df_absplice_rna = pd.read_csv(snakemake.input['absplice_rna'])

index = ['variant', 'gene_id', 'tissue']
rna_cols = ['sample', 'delta_psi_cat', 'pValueGene_g_minus_log10', 'AbSplice_RNA']
df_absplice = df_absplice_rna.set_index(index)[rna_cols].join(
    df_absplice_dna.set_index(index), 
    how='outer',
).reset_index()

cols = [
    'variant', 'gene_id', 'tissue', 'sample',
    'AbSplice_DNA', 'AbSplice_RNA',
    'delta_psi_cat', 'pValueGene_g_minus_log10', 
    'delta_logit_psi', 'delta_psi', 'delta_score',
    'splice_site_is_expressed',   
    'ref_psi', 'median_n',
    'junction', 'splice_site', 'event_type',
    'acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss', 
    'acceptor_gain_position', 'acceptor_loss_position', 'donor_gain_position', 'donor_loss_position'
]
df_absplice[cols].to_csv(snakemake.output['absplice'], index=False)