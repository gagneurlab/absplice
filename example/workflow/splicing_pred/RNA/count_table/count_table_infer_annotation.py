import pandas as pd
from splicemap import SpliceCountTable as CountTable

ct = CountTable.read_csv(snakemake.input['count_table'], 
                        name=snakemake.params['tissue_cat'])

gtf_file = snakemake.input['gtf_file']
ct.infer_annotation(gtf_file)

df_annotation = ct.annotation
df_annotation = df_annotation[
    df_annotation['gene_type'].str.contains('protein_coding')
]
df_annotation = df_annotation.drop(columns=['transcript_id', 'gene_name', 'gene_type']).reset_index()
df_annotation['gene_id'] = df_annotation['gene_id'].apply(lambda x: x.split(';'))
df_annotation = df_annotation.explode('gene_id')

coding_genes = pd.read_csv(snakemake.input['coding_genes'])
df_annotation = df_annotation.set_index('gene_id').join(
    coding_genes.set_index('gene_id')[['gene_name', 'gene_type']]
    ).reset_index()
df_annotation[
    df_annotation['gene_type'] == 'protein_coding'
]

if 'chr' not in df_annotation['junctions'].values[0]:
    df_annotation['junctions'] = df_annotation['junctions'].apply(lambda x: 'chr' + x)

df_annotation.to_csv(snakemake.output['count_table_with_annotation'], index=False)
