from absplice import SplicingOutlierResult

splicing_result = SplicingOutlierResult(
        df_mmsplice=snakemake.input['df_mmsplice'], 
        df_spliceai=snakemake.input['df_spliceai'], 
        df_mmsplice_cat=snakemake.input['df_mmsplice_cat'], 
        df_outliers_cat=snakemake.input['df_outliers_cat'], 
        df_var_samples=snakemake.input['var_samples_df'], 
    )
splicing_result.predict_absplice_rna()
df_absplice_rna = splicing_result._absplice_rna

df_absplice_rna = df_absplice_rna.rename(columns={
    'pValueGene_g_minus_log10': 'FRASER_pval_cat'
})

df_absplice_rna.to_csv(snakemake.output['absplice_rna'])