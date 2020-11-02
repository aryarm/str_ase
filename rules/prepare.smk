import csv
import warnings
from pathlib import Path
from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.26.1")

configfile: "rules/prepare.yml"
configfile: "config.yml"



def check_config(value, default=False, place=config):
    """ return true if config value exists and is true """
    return place[value] if (value in place and place[value]) else default

# set the output directory if it isn't set already
config['data_dir'] = check_config('data', default='data')
# get the current chr
chr_name = config['region_hg19'].split(':')[0]

def read_samples():
    """
        Function to get names and dna fastq paths from a sample file
        specified in the configuration. Input file is expected to have 3
        columns: <unique_sample_id> <fastq1_path> <fastq2_path> or
        <unique_sample_id> <paired_bam_path> <bed_path>. Modify this function
        as needed to provide a dictionary of sample_id keys and either a tuple
        of strings: (fastq1, fastq2) OR a single string: paired_bam
    """
    samp_dict = {}
    for line in csv.reader(open(config['samples_file']), delimiter="\t"):
        if len(line) == 3:
            samp_dict[line[1]] = [line[0], line[2]]
        else:
            raise ValueError('Your samples_file is not formatted correctly. Make sure that it has the correct number of tab-separated columns for every row.')
    return samp_dict
SAMP = read_samples()

# the user can change config['SAMP_NAMES'] here (or define it in the config
# file) to contain whichever sample names they'd like to run the pipeline on
if 'SAMP_NAMES' not in config or not config['SAMP_NAMES']:
    config['SAMP_NAMES'] = list(SAMP.keys())
else:
    # double check that the user isn't asking for samples they haven't provided
    user_samps = set(config['SAMP_NAMES'])
    config['SAMP_NAMES'] = list(set(SAMP.keys()).intersection(user_samps))
    if len(config['SAMP_NAMES']) != len(user_samps):
        warnings.warn("Not all of the samples requested have provided input. Proceeding with as many samples as is possible...")



rule all:
    input:
        config['str_vcf'], config['snp_vcf'],
        config['ref_genome'], config['ref_panel'],
        config['lift_over_chain']

rule create_ref_genome:
    input:
        ref_genome = config['ref_genome_hg19']
    params:
        chr_name = chr_name
    output:
        ref_genome = config['ref_genome'],
        ref_genome_idx = config['ref_genome']+".fai"
    conda: "../envs/htslib.yml"
    shell:
        "samtools faidx {input.ref_genome} {params.chr_name} > {output.ref_genome} && "
        "samtools faidx {output.ref_genome}"

# TODO: make this download from ucsc if the input doesn't exist
rule create_lift_over:
    input:
        chain = config['lift_over']
    output:
        chain = config['lift_over_chain']
    conda: "../envs/default.yml"
    shell:
        "ln -sf {input.chain} {output.chain}"

# TODO: detect whether index exists automatically and skip this step
rule index_snp_source_vcf:
    input:
        snp_source_vcf = config['snp_vcfs']
    output:
        link = temp(config['data_dir']+"/snp-full.vcf.gz"),
        index = temp(config['data_dir']+"/snp-full.vcf.gz.tbi")
    conda: "../envs/htslib.yml"
    shell:
        "ln -sf {input.snp_source_vcf} {output.link} && "
        "tabix -p vcf {output.link}"

rule create_snp_vcf:
    input:
        source = rules.index_snp_source_vcf.output.link,
        source_idx = rules.index_snp_source_vcf.output.index,
        samples = config['samples_file']
    params:
        region = "chr"+config['region_hg38']
    output:
        vcf = config['snp_vcf'],
        vcf_idx = config['snp_vcf']+".tbi"
    conda: "../envs/htslib.yml"
    shell:
        "bcftools view -r '{params.region}' -S <(cut -f1 {input.samples}) {input.source} -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule reheader_original_str_vcf:
    input:
        vcf = config['str_vcfs'],
        ref_genome_idx = config['ref_genome_hg19']+".fai"
    params:
        new_samp_name = lambda wildcards: SAMP[wildcards.sample][0]
    output:
        vcf = temp(config['data_dir']+"/reheader_strs/{sample}.vcf.gz"),
        vcf_idx = temp(config['data_dir']+"/reheader_strs/{sample}.vcf.gz.tbi")
    conda: "../envs/htslib.yml"
    shell:
        "bcftools reheader -f {input.ref_genome_idx} -s <(echo '{params.new_samp_name}') {input.vcf} -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule merge_full_strs:
    input:
        vcfs = expand(rules.reheader_original_str_vcf.output.vcf, sample=config['SAMP_NAMES'])
    params:
        vcfs = lambda wildcards, input: ",".join(input.vcfs),
        unzipped_vcf = lambda wildcards, output: Path(output.vcf).with_suffix(''),
        vcf_prefix = lambda wildcards, output: Path(Path(output.vcf).with_suffix('')).with_suffix('')
    output:
        vcf = temp(config['data_dir']+"/str-full.vcf.gz")
    conda: "../envs/default.yml"
    shell:
        "mergeSTR --vcftype hipstr --vcfs '{params.vcfs}' --out {params.vcf_prefix} && "
        "bgzip {params.unzipped_vcf}"

rule reheader_full_str_vcf:
    input:
        ref_genome_idx = config['ref_genome_hg19']+".fai",
        full_vcf = rules.merge_full_strs.output.vcf
    output:
        vcf = temp(config['data_dir']+"/str-full.reheader.vcf.gz")
    conda: "../envs/htslib.yml"
    shell:
        "bcftools reheader -f {input.ref_genome_idx} {input.full_vcf} -o {output.vcf}"

rule sort_full_str_vcf:
    input:
        reheader_vcf = rules.reheader_full_str_vcf.output.vcf
    params:
        data_dir = config['data_dir']
    output:
        vcf = temp(config['data_dir']+"/str-full.sorted.vcf.gz"),
        vcf_idx = temp(config['data_dir']+"/str-full.sorted.vcf.gz.tbi")
    conda: "../envs/htslib.yml"
    shell:
        "bcftools sort -Oz -o {output.vcf} {input.reheader_vcf} -T {params.data_dir} && "
        "tabix -p vcf {output.vcf}"

rule create_str_vcf:
    input:
        sorted_vcf = rules.sort_full_str_vcf.output.vcf,
        sorted_vcf_idx = rules.sort_full_str_vcf.output.vcf_idx
    params:
        region = config['region_hg19']
    output:
        vcf = config['str_vcf'],
        vcf_idx = config['str_vcf']+".tbi"
    conda: "../envs/htslib.yml"
    shell:
        "bcftools view -r '{params.region}' {input.sorted_vcf} -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule create_ref_panel:
    input:
        ref_panel = config['ref_panel_dir']
    output:
        ref_panel = directory(config['ref_panel'])
    conda: "../envs/default.yml"
    shell:
        "ln -sfn {input.ref_panel} {output.ref_panel}; "
        "test -L {output.ref_panel} && test -d {output.ref_panel}"