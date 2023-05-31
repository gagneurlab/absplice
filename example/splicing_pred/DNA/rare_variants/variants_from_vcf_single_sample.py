import pandas as pd
from kipoiseq.extractors.vcf_seq import MultiSampleVCF
from gnomad_rocksdb import GnomadMafDB
from absplice.utils import VariantMafFilter, PrivateVariantFilter, ReadDepthFilter, GQFilter

maf = GnomadMafDB(snakemake.input['maf'])
vcf = MultiSampleVCF(snakemake.input['vcf'])

# filter vcf based on provided filters
vcf_filtered = vcf.query_all()

vcf_filtered = vcf_filtered.filter(lambda v: v.filter is None)

if snakemake.params['filter_maf']:
    vcf_filtered = vcf_filtered.filter(
       lambda v: VariantMafFilter(cutoff=float(snakemake.params['maf_cutoff']), population=maf)(v))

# if snakemake.params['filter_genotype_qual']:
#     vcf_filtered = vcf_filtered.filter(
#        lambda v: GQFilter(vcf, min_GQ=int(snakemake.params['genotype_quality']))(v))
    
# write variants into csv file    
variants = [x.__str__() for x in vcf_filtered]

# DNA ID
with open(snakemake.input['vcf_DNA_ID']) as f:
    lines = f.readlines()
assert len(lines) == 1
DNA_ID = lines[0]

df = pd.DataFrame({
    'variant': variants,
    'sample': DNA_ID,
})

df.to_csv(snakemake.output['filtered_variants'], index=False)
