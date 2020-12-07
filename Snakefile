import csv
from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.27.3")

configfile: "config.yml"



def check_config(value, default=False, place=config):
    """return true if config value exists and is true"""
    return place[value] if (value in place and place[value]) else default

# set the output directory if it isn't set already
config['out'] = check_config('out', default='out')
# and set the SAMP_NAMES if they haven't been set yet
config['SAMP_NAMES'] = check_config('SAMP_NAMES', default=[])



rule all:
    input:
        config['out']+"/str_counts.tsv.gz"

rule lift_over:
    """lift SNP VCF from hg38 to hg19 using Picard's LiftoverVcf"""
    input:
        vcf = config['snp_vcf'],
        idx = config['snp_vcf']+".tbi",
        chain = config['lift_over_chain'],
        ref = config['ref_genome_hg19']
    output:
        vcf = temp(config['out']+"/snp.hg19.vcf.gz"),
        vcf_unmapped = temp(config['out']+"/snp.hg19.unmap.vcf.gz")
    conda: "envs/default.yml"
    shell:
        "gatk LiftoverVcf -I {input.vcf} -O {output.vcf} -C {input.chain} --REJECT {output.vcf_unmapped} -R {input.ref} -WMC true -LFI false"

rule remove_dups:
    """remove duplicate variants in the SNP VCF from the liftover"""
    input:
        vcf = rules.lift_over.output.vcf,
        ref = config['ref_genome_hg19']
    output:
        vcf = temp(config['out']+"/snp.hg19.rmdup.vcf.gz"),
        idx = temp(config['out']+"/snp.hg19.rmdup.vcf.gz.tbi")
    conda: "envs/htslib.yml"
    shell:
        "bcftools norm -Oz -o {output.vcf} -d all -cx -f {input.ref} {input.vcf} && "
        "tabix -p vcf {output.vcf}"

rule hg192b37:
    """remove chr prefixes from SNP VCF"""
    input:
        vcf = rules.remove_dups.output.vcf
    output:
        vcf = config['out']+"/snp.vcf.gz",
        idx = config['out']+"/snp.vcf.gz.tbi"
    conda: "envs/htslib.yml"
    shell:
        "zcat {input.vcf} | "
        "sed 's/^chr//; s/^##contig=<ID=chr/##contig=<ID=/' | "
        "bgzip > {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule vcf_merge:
    input:
        str_vcf = config['str_vcf'],
        snp_vcf = rules.hg192b37.output.vcf
    output:
        vcf = temp(config['out']+"/merged.vcf.gz"),
        idx = temp(config['out']+"/merged.vcf.gz.tbi")
    conda: "envs/htslib.yml"
    shell:
        "bcftools concat -Ob {input} | "
        "bcftools sort -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule remove_empty_alts:
    input:
        vcf = rules.vcf_merge.output.vcf,
        idx = rules.vcf_merge.output.idx
    output:
        vcf = temp(config['out']+"/merged.noalts.vcf.gz"),
        idx = temp(config['out']+"/merged.noalts.vcf.gz.tbi")
    conda: "envs/htslib.yml"
    shell:
        "scripts/remove_empty_alts.py {input.vcf} | bgzip > {output.vcf} && "
        "tabix -p vcf {output.vcf}"

checkpoint vcf_chroms:
    """get the chroms from a VCF"""
    input:
        vcf = rules.remove_empty_alts.output.vcf,
        vcf_index = rules.remove_empty_alts.output.idx
    output:
        chroms = config['out'] + "/merged/chroms.txt"
    conda: "envs/htslib.yml"
    shell:
        "tabix --list-chroms {input.vcf} > {output.chroms}"

rule split_vcf_by_chr:
    """Split the provided VCF file by chromosome and bgzip it"""
    input:
        vcf = rules.remove_empty_alts.output.vcf,
        vcf_index = rules.remove_empty_alts.output.idx,
    output:
        vcf = config['out'] + "/merged/{chr}.vcf.gz",
        idx = config['out'] + "/merged/{chr}.vcf.gz.tbi"
    conda: "envs/htslib.yml"
    shell:
        "tabix -h {input.vcf} {wildcards.chr} | bgzip > {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule conform_gt:
    input:
        vcf = rules.split_vcf_by_chr.output.vcf,
        ref =  config['ref_panel']+"/1kg.snp.str.chr{chr}.vcf.gz",
        conform_gt = config['conform_gt_jar']
    params:
        region = lambda wildcards: wildcards.chr,
        vcf_prefix = lambda wildcards, output: output.vcf[:-len(".vcf.gz")]
    output:
        vcf = config['out']+"/consistent/snp.str.chr{chr}.vcf.gz",
        log = temp(config['out']+"/consistent/snp.str.chr{chr}.log")
    conda: "envs/default.yml"
    shell:
        "java -jar {input.conform_gt} gt={input.vcf} ref={input.ref} chrom={params.region} match=POS out={params.vcf_prefix}"

rule beagle:
    input:
        gt = rules.conform_gt.output.vcf,
        ref = config['ref_panel']+"/1kg.snp.str.chr{chr}.vcf.gz",
        genetic_map = config['genetic_map']+"/plink.chr{chr}.GRCh37.map",
        beagle = config['beagle_jar']
    params:
        region = lambda wildcards: wildcards.chr,
        vcf_prefix = lambda wildcards, output: output.vcf[:-len(".vcf.gz")]
    output:
        vcf = temp(config['out']+"/phased/snp.str.chr{chr}.ungz.vcf.gz"),
        log = temp(config['out']+"/phased/snp.str.chr{chr}.ungz.log")
    conda: "envs/default.yml"
    shell:
        "java -jar {input.beagle} gt={input.gt} ref={input.ref} out={params.vcf_prefix} map={input.genetic_map} chrom={params.region}"

rule index_beagle:
    """properly bgzip the output of beagle and add contigs to the header"""
    input:
        vcf = rules.beagle.output.vcf,
        ref = config['ref_genome']+".fai"
    output:
        vcf = config['out']+"/phased/snp.str.chr{chr}.vcf.gz",
        idx = config['out']+"/phased/snp.str.chr{chr}.vcf.gz.tbi",
    conda: "envs/htslib.yml"
    shell:
        "bcftools reheader -o {output.vcf} -f {input.ref} {input.vcf} && "
        "tabix -p vcf {output.vcf}"

def get_split_vcf(wildcards):
    """get the output of index_beagle expanded for every chromosome"""
    with checkpoints.vcf_chroms.get().output.chroms.open() as f:
        return expand(
            rules.index_beagle.output.vcf,
            chr=filter(lambda x: len(x), f.read().split('\n'))
        )

rule merge_beagle:
    """merge the beagle output across all chromosomes"""
    input:
        vcfs = get_split_vcf
    output:
        vcf = config['out']+"/phased.vcf.gz",
        idx = config['out']+"/phased.vcf.gz.tbi"
    conda: "envs/htslib.yml"
    shell:
        "bcftools concat -Ob -l {input} | "
        "bcftools sort -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule lift_counts:
    input:
        counts = config['snp_counts'],
        chain = config['lift_over_chain']
    output:
        counts = config['out']+"/snp_counts/{sample}.tsv.gz"
    conda: "envs/default.yml"
    shell:
        "scripts/liftover_counts.py -i {input.chain} {input.counts} | "
        "sed 's/^chr//' | "
        "gzip > {output.counts}"

rule str_counts:
    """get per sample read counts of STRs from surrounding SNPs"""
    input:
        vcf = rules.merge_beagle.output.vcf,
        idx = rules.merge_beagle.output.idx,
        counts = rules.lift_counts.output.counts
    output:
        counts = config['out']+"/str_counts/{sample}.tsv.gz"
    conda: "envs/htslib.yml"
    shell:
        "scripts/as_counts.py {input.counts} {wildcards.sample} "
        "{input.vcf} | gzip > {output.counts}"

checkpoint vcf_samples:
    """get the sample IDs from the merged VCF"""
    input:
        vcf = rules.merge_beagle.output.vcf,
        vcf_index = rules.merge_beagle.output.idx
    output:
        samples = config['out'] + "/phased/samples.txt"
    conda: "envs/htslib.yml"
    shell:
        "bcftools query -l {input.vcf} > {output.samples}"

def get_split_counts(wildcards):
    """get the output of str_counts for every sample"""
    with checkpoints.vcf_samples.get().output.samples.open() as f:
        return expand(
            rules.str_counts.output.counts,
            sample=filter(lambda x: len(x), f.read().split('\n'))
        )

rule merge_str_counts:
    """merge str_counts TSVs across all samples and sort on str, tissue, subject"""
    input:
        counts = get_split_counts
    output:
        merged = config['out']+"/str_counts.tsv.gz"
    conda: "envs/default.yml"
    shell:
        "zcat {input.counts} | ("
        "read -r head && echo \"$head\"; "
        "grep -Fxv \"$head\""
        ") | gzip > {output.merged}"

rule prioritize:
    """compute an importance metric for each STR and then sort by that metric"""
    input:
        counts = rules.merge_str_counts.output.merged
    output:
        counts = config['out']+"/str_counts.sort.tsv.gz"
    conda: "envs/htslib.yml"
    shell:
        "scripts/prioritize_STRs.py {input} {output}"
