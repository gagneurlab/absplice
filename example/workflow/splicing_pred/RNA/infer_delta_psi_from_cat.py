import pdb
import pandas as pd
from absplice.result import SplicingOutlierResult
from absplice.cat_dataloader import CatInference
import os

# if SpliceMaps for accessible tissue are provided reference PSI estimates are taken from them
# (this can give a better estimate, especially if the cohort of provided RNA-seq samples is small)
if os.path.exists(snakemake.input['splicemap_cat5']):
    splicemap_cat5 = snakemake.input['splicemap_cat5']
else:
    splicemap_cat5 = None
if os.path.exists(snakemake.input['splicemap_cat3']):
    splicemap_cat3 = snakemake.input['splicemap_cat3']
else:
    splicemap_cat3 = None

# Infer CAT
cat_dl = CatInference(
    splicemap5=list(snakemake.input['splicemap_5']),
    splicemap3=list(snakemake.input['splicemap_3']),
    splicemap_cat5=splicemap_cat5,
    splicemap_cat3=splicemap_cat3,
    count_cat=snakemake.input['cat_count_table'],
    name=snakemake.params['tissue_cat']
)

result = SplicingOutlierResult(
    df_mmsplice = snakemake.input['mmsplice_splicemap'],
)
result.add_samples(pd.read_csv(snakemake.input['var_samples_df']))

result.infer_cat(cat_dl, progress=True)
result.df_mmsplice_cat.to_csv(snakemake.output['result'])

