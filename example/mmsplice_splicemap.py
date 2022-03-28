from splicing_outlier_prediction import SpliceOutlier, SpliceOutlierDataloader

dl = SpliceOutlierDataloader(
    snakemake.input['fasta'], snakemake.input['vcf'],
    splicemap5=list(snakemake.params['splicemap_5']),
    splicemap3=list(snakemake.params['splicemap_3'])
)

model = SpliceOutlier()
model.predict_save(dl, snakemake.output['result'])