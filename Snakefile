from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.26.1")

configfile: "config.yml"



def check_config(value, default=False, place=config):
    """ return true if config value exists and is true """
    return place[value] if (value in place and place[value]) else default

# set the output directory if it isn't set already
config['out'] = check_config('out', default='out')
# and set the SAMP_NAMES if they haven't been set yet
config['SAMP_NAMES'] = check_config('SAMP_NAMES', default=[])



rule all:
	input: config['out']+"/phased/snp.str.chr1.vcf.gz"

rule beagle:
	input:
		# TODO: add a rule for splitting the VCF by chrom + use the chrom-specific file
		gt = config['vcf'],
		ref = config['ref_panel']+"/1kg.snp.str.chr{chr}.vcf.gz"
	output:
		phased = config['out']+"/phased/snp.str.chr{chr}.vcf.gz"
	conda: "env.yml"
	shell:
		"beagle gt={input.gt} ref={input.ref} out={output.phased}"
