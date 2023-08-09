from absplice import SplicingOutlierResult

splicing_result = SplicingOutlierResult(
        df_mmsplice=snakemake.input['mmsplice_splicemap'], 
        df_spliceai=snakemake.input['spliceai'],
    )
splicing_result.predict_absplice_dna()
splicing_result._absplice_dna.to_csv(snakemake.output['absplice_dna'])