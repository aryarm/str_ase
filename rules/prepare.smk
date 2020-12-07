import os
import csv
import warnings
from pathlib import Path
from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.27.3")

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
    samp_dict = {}
    for line in csv.reader(open(config['samples_file']), delimiter="\t"):
        if len(line) == 2:
            samp_dict[line[1]] = line[0]
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
        config['ref_genome'], os.path.splitext(config['ref_genome'])[0]+".dict",
        config['ref_panel']+"/1kg.snp.str.chr"+chr_name+".vcf.gz",
        config['lift_over_chain'], expand(
            config['snp_counts'],
            sample=[SAMP[samp] for samp in config['SAMP_NAMES']]
        ),
        expand(config['genetic_map']+"/plink.chr{chr}.GRCh37.map", chr=chr_name)

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

rule create_ref_dict:
    input:
        ref_genome = config['ref_genome']
    output:
        ref_genome_dict = os.path.splitext(config['ref_genome'])[0]+".dict"
    conda: "../envs/default.yml"
    shell:
        "gatk CreateSequenceDictionary -R {input.ref_genome}"

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
        samples = config['samples_file'],
        ref = config['ref_genome_hg38']
    params:
        region = "chr"+config['region_hg38']
    output:
        vcf = config['snp_vcf'],
        vcf_idx = config['snp_vcf']+".tbi"
    conda: "../envs/default.yml"
    shell:
        "bcftools view -r '{params.region}' -S <(cut -f1 {input.samples}) {input.source} -Ov -o- | "
        #"vcf-convert -v 4.1 -r {input.ref} | " # TODO: find some way of converting from vcf v4.2 to v4.1 since this doesn't work
        "grep -Pv '^##contig=<ID=.*(?<!(chr2)),length=' | " # remove all other contig lines from the header
        "bgzip > {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule reheader_original_str_vcf:
    input:
        vcf = config['str_vcfs'],
        ref_genome_idx = config['ref_genome_hg19']+".fai"
    params:
        new_samp_name = lambda wildcards: SAMP[wildcards.sample]
    output:
        vcf = temp(config['data_dir']+"/reheader_strs/{sample}.vcf.gz"),
        vcf_idx = temp(config['data_dir']+"/reheader_strs/{sample}.vcf.gz.tbi")
    conda: "../envs/htslib.yml"
    shell:
        "bcftools reheader -f {input.ref_genome_idx} -s <(echo '{params.new_samp_name}') {input.vcf} -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule subset_original_str_vcf:
    input:
        str_vcf = rules.reheader_original_str_vcf.output.vcf,
        str_vcf_idx = rules.reheader_original_str_vcf.output.vcf_idx
    params:
        region = config['region_hg19']
    output:
        vcf = temp(config['data_dir']+"/subset_strs/{sample}.vcf.gz"),
        vcf_idx = temp(config['data_dir']+"/subset_strs/{sample}.vcf.gz.tbi")
    conda: "../envs/htslib.yml"
    shell:
        "bcftools view -r '{params.region}' {input.str_vcf} -Oz -o {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule merge_full_strs:
    input:
        vcfs = expand(rules.subset_original_str_vcf.output.vcf, sample=config['SAMP_NAMES']),
        vcf_idxs = expand(rules.subset_original_str_vcf.output.vcf_idx, sample=config['SAMP_NAMES'])
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
        "bcftools sort -Oz -o {output.vcf} {input.reheader_vcf} -T {params.data_dir}"
        " && tabix -p vcf {output.vcf}"

rule unphase_full_str_vcf:
    input:
        vcf = rules.sort_full_str_vcf.output.vcf,
        vcf_idx = rules.sort_full_str_vcf.output.vcf_idx,
        ref_idx = rules.create_ref_genome.output.ref_genome_idx
    output:
        vcf = temp(config['data_dir']+"/str-full.unphased.vcf.gz"),
        vcf_idx = temp(config['data_dir']+"/str-full.unphased.vcf.gz.tbi")
    conda: "../envs/htslib.yml"
    shell:
        "scripts/unphase.py {input.vcf} | bgzip > {output.vcf} && "
        "tabix -p vcf {output.vcf}"

rule reheader_final_str_vcf:
    input:
        vcf = rules.unphase_full_str_vcf.output.vcf,
        ref_idx = rules.create_ref_genome.output.ref_genome_idx
    output:
        vcf = config['str_vcf'],
        vcf_idx = config['str_vcf']+".tbi"
    conda: "../envs/htslib.yml"
    shell:
        "bcftools reheader -o {output.vcf} -f {input.ref_idx} {input.vcf} && "
        "tabix -p vcf {output.vcf}"

rule create_ref_panel:
    input:
        ref_panel = config['ref_panel_dir']+"/1kg.snp.str.chr"+chr_name+".vcf.gz",
        ref_panel_idx = config['ref_panel_dir']+"1kg.snp.str.chr"+chr_name+".vcf.gz.tbi"
    params:
        region = config['region_hg19']
    output:
        ref_panel = config['ref_panel']+"/1kg.snp.str.chr"+chr_name+".vcf.gz",
        ref_panel_idx = config['ref_panel']+"/1kg.snp.str.chr"+chr_name+".vcf.gz.tbi"
    conda: "../envs/htslib.yml"
    shell:
        "bcftools view -r '{params.region}' {input.ref_panel} -Oz -o {output.ref_panel}"
        " && tabix -p vcf {output.ref_panel}"

rule download_plink_map:
    params:
        url = "http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip",
        out_dir = lambda wildcards, output: Path(output[0]).parent
    output:
        temp(config['genetic_map']+"/plink.GRCh37.map.zip"),
        expand(
            config['genetic_map']+"/plink.chr{chr}.GRCh37.map",
            chr=chr_name
        )
    conda: "../envs/default.yml"
    shadow: "minimal"
    shell:
        "wget -P {params.out_dir} {params.url} && "
        "unzip -d {params.out_dir} {output[0]}"

rule create_snp_counts:
    input:
        snp_counts = config['snp_counts_full']
    params:
        chr_name = "chr"+chr_name,
        start = config['region_hg38'].split(':')[1].split('-')[0],
        end = config['region_hg38'].split(':')[1].split('-')[1],
    output:
        snp_counts = config['snp_counts']
    conda: "../envs/default.yml"
    shell:
        "zcat {input.snp_counts} | ("
        "read -r head && echo \"$head\"; "
        "awk -F '\\t' -v 'OFS=\\t' "
        "'$1 == \"{params.chr_name}\" && $2 < {params.end} && $2 > {params.start}'"
        ") | gzip > {output.snp_counts}"
