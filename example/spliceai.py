from spliceai_pipeline.spliceAI import SpliceAI

genome_map = {
    'hg19': 'grch37',
    'hg38': 'grch38'
}

if snakemake.params['lookup_only']:
    model = SpliceAI(db_path=snakemake.input['db'])
else:
    model = SpliceAI(fasta=snakemake.input['fasta'],
                     annotation=genome_map[snakemake.params['genome']],
                     db_path=snakemake.input['db'])

model.predict_save(snakemake.input['vcf'],
                   snakemake.output['result'])