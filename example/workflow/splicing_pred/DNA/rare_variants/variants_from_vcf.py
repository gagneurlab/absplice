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

if snakemake.params['filter_private']:
    vcf_filtered = vcf_filtered.filter(
       lambda v: PrivateVariantFilter(vcf, max_num_samples=int(snakemake.params['max_num_samples']))(v))
    
# if snakemake.params['filter_genotype_qual']:
#     vcf_filtered = vcf_filtered.filter(
#        lambda v: GQFilter(vcf, min_GQ=int(snakemake.params['genotype_quality']))(v))
    
vcf_filtered.to_sample_csv(snakemake.output['filtered_variants'])
