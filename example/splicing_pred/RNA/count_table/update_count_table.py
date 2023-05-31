import pandas as pd
from splicemap import SpliceCountTable as CountTable

# Read in sample mapping from RNA_ID to DNA_ID
df = pd.read_csv(snakemake.input['DROP_annotation'], sep='\t')
# subset for accessible tissue
df = df[
    df['DROP_GROUP'].apply(lambda x: x.split(','))\
        .apply(lambda x: snakemake.wildcards['tissue_cat'] in x)
]
sampleID_map_RNA_to_DNA = pd.Series(df['DNA_ID'].values,index=df['sampleID']).to_dict()

# Read and convert to pyranges
ct = CountTable.read_csv(snakemake.input['raw_count_table'],
                         name=snakemake.wildcards['tissue_cat'])

fasta_file = snakemake.input['fasta_file']
ct.infer_strand(fasta_file, progress=True)

ct.update_samples(sampleID_map_RNA_to_DNA)

if 'chr' not in ct.df['Chromosome'].values[0]:
    ct.df['Chromosome'] = 'chr' + ct.df['Chromosome'].astype(str)

ct.to_csv(snakemake.output['updated_count_table'])
