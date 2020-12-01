#!/usr/bin/env python
import sys
import argparse
import pandas as pd
from pysam import VariantFile



def get_hets(vcf, sample):
    """extract heterozygous STRs and SNPs from the vcf"""
    vcf = VariantFile(vcf)
    vcf.subset_samples([sample])
    for rec in vcf.fetch():
        variant = rec.samples[sample]
        is_het = variant['GT'][0] != variant['GT'][1]
        is_snp = not ((len(variant.alleles[0]) - len(variant.alleles[1])))
        # limit our analysis to only variants that are heterozygous
        # and they must be either STRs or SNPs
        if (is_het and (is_snp or rec.id.startswith("STR_"))):
            yield rec


def get_str_snp_pairs(variants):
    """return tuples of (str, nearest snp)"""
    last_snp = None
    # iterate through the SNPs
    while True:
        strs = []
        # iterate through the strs
        for variant in variants:
            if not variant.id.startswith("STR_"):
                break
            strs.append(variant)
        # check is the final variant an STR?
        # if it is the last het variant in the VCF, it might be an STR
        if variant.id.startswith("STR_"):
            for STR in strs:
                yield (STR, last_snp)
            return
        # check: did we see any strs?
        if len(strs):
            # Now, we must find the nearest SNP
            # If we didn't previously see a SNP,
            # this is just going to be the current variant
            if last_snp is None:
                for STR in strs:
                    yield (STR, variant)
            else:
                # we must decide whether to use the current SNP or the last one
                snps = (last_snp, variant)
                for STR in strs:
                    # TODO: also consider CHROM, not just POS
                    distances = (STR.pos-last_snp.pos, variant.pos-STR.pos)
                    yield (STR, snps[distances.index(min(distances))])
        else:
            # just save the current SNP and move on
            last_snp = variant


def get_snp_counts(snp_counts):
    """parse the SNP ASE table"""
    return pd.read_csv(snp_counts, sep="\t", header=0, index_col=['CHR', 'POS'])


def get_genes(genes):
    """parse the gene annotations file"""
    pass


def get_str_counts(pairs, snp_counts, genes, sample):
    for pair in pairs:
        STR, snp = pair
        if (snp.chrom, snp.pos) in snp_counts.index:
            yield (STR, snp, snp_counts[(snp.chrom, snp.pos)])
        # discard the current STR
        # TODO: try to find a new SNP?


def write_out(out, sample, str_counts):
    print("ahhh this should be a TSV with appropriately named columns", file=out)
    # STR, gene, sample, ASE_SNP,_location, STR a1|a2, SNP a1|a2, SNP ASE stats


def get_args():
    """parse the arguments to this script using argparse"""
    parser = argparse.ArgumentParser(description=
        "Get per-sample counts of reads overlapping STRs "
        "by using counts from nearby SNPs."
    )
    parser.add_argument(
        '-o', '--out', help='output file name (in tsv format)',
        type=argparse.FileType('w')
    )
    parser.add_argument(
        'vcf', help="sorted and phased STRs (prefixed with 'STR_') and SNPs"
    )
    parser.add_argument(
        'counts', help='WASP-corrected counts of SNPs for this sample'
        # TODO: describe the structure of this tsv file better
    )
    parser.add_argument(
        'genes', help='GENCODE gene annotations'
        # TODO: describe the structure of this tsv file better
    )
    parser.add_argument(
        'sample', help='the name of the sample for which to extract STR counts'
    )
    args = parser.parse_args()
    return args


def main(args):
    write_out(
        args.out, args.sample,
        get_str_counts(
            get_str_snp_pairs(get_hets(args.vcf, args.sample)),
            get_snp_counts(args.counts),
            get_genes(args.genes),
            args.sample
        )
    )


if __name__ == "__main__":
    main(get_args())
