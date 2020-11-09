from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.27.3")

configfile: "config.yml"



def check_config(value, default=False, place=config):
    """ return true if config value exists and is true """
    return place[value] if (value in place and place[value]) else default

# set the output directory if it isn't set already
config['out'] = check_config('out', default='out')
# and set the SAMP_NAMES if they haven't been set yet
config['SAMP_NAMES'] = check_config('SAMP_NAMES', default=[])



rule all:
    input: config['out']+"/phased/snp.str.chr2.vcf.gz"

rule lift_over:
    """ lift SNP VCF from hg38 to hg19 using CrossMap """
    input:
        vcf = config['snp_vcf'],
        chain = config['lift_over_chain'],
        ref = config['ref_genome']
    output:
        vcf = temp(config['out']+"/snp.hg19.vcf"),
        vcf_unmapped = temp(config['out']+"/snp.hg19.vcf.unmap")
    conda: "envs/crossmap.yml"
    shell:
        "CrossMap.py vcf {input.chain} {input.vcf} {input.ref} {output.vcf}"

rule remove_chr_prefix:
    """ remove chr prefixes from SNP VCF """
    input:
        vcf = rules.lift_over.output.vcf
    output:
        vcf = config['out']+"/snp.hg19.vcf.gz"
    conda: "envs/default.yml"
    shell:
        "sed 's/^chr//; s/^##contig=<ID=chr/##contig=<ID=/' {input.vcf} | "
        "bgzip > {output.vcf}"

rule vcf_merge:
    input:
        str_vcf = config['str_vcf'],
        snp_vcf = rules.remove_chr_prefix.output.vcf
    output:
        vcf = temp(config['out']+"/merged.vcf.gz"),
        idx = temp(config['out']+"/merged.vcf.gz.tbi")
    conda: "envs/htslib.yml"
    shell:
        "bcftools concat -Ob {input} | "
        "bcftools sort -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"

checkpoint vcf_chroms:
    """get the chroms from a VCF"""
    input:
        vcf = rules.vcf_merge.output.vcf,
        vcf_index = rules.vcf_merge.output.idx
    output:
        chroms = config['out'] + "/merged/chroms.txt"
    conda: "envs/default.yml"
    shell:
        "tabix --list-chroms {input.vcf} > {output.chroms}"

rule split_vcf_by_chr:
    """Split the provided VCF file by chromosome and bgzip it"""
    input:
        vcf = rules.vcf_merge.output.vcf,
        vcf_index = rules.vcf_merge.output.idx,
    output:
        vcf = config['out'] + "/merged/{chr}.vcf.gz",
        idx = config['out'] + "/merged/{chr}.vcf.gz.tbi"
    conda: "envs/default.yml"
    shell:
        "tabix -h {input.vcf} {wildcards.chr} | bgzip > {output.vcf} && "
        "tabix -p vcf {output.vcf}"

def get_split_vcf(wildcards):
    with checkpoints.vcf_chroms.get().output.chroms.open() as f:
        return expand(
            rules.split_vcf_by_chr.output,
            chr=filter(lambda x: len(x), f.read().split('\n'))
        )

rule conform_gt:
    input:
        vcf = rules.split_vcf_by_chr.output.vcf,
        ref =  config['ref_panel']+"/1kg.snp.str.chr{chr}.vcf.gz"
    params:
        region = lambda wildcards: wildcards.chr
    output:
        vcf = config['out']+"/consistent/snp.str.chr{chr}.vcf.gz"
    conda: "envs/default.yml"
    shell:
        "conform-gt gt={input.vcf} ref={input.ref} chrom={params.region} match=POS out={output.vcf}"

rule beagle:
    input:
        gt = rules.conform_gt.output.vcf,
        ref = config['ref_panel']+"/1kg.snp.str.chr{chr}.vcf.gz"
    output:
        vcf = config['out']+"/phased/snp.str.chr{chr}.vcf.gz"
    conda: "envs/default.yml"
    shell:
        "beagle gt={input.gt} ref={input.ref} out={output.vcf}"
