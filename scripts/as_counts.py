#!/usr/bin/env python
import sys
import argparse
from pysam import VariantFile



def get_hets(vcf, sample):
    """parse heterozygous variants from the vcf"""
    vcf = VariantFile(vcf)
    vcf.subset_samples([sample])
    for rec in vcf.fetch():
        # limit our analysis to only variants that are heterozygous
        # TODO: are hets just 0|1 or 1|0? or can they also be defined like this?
        if rec.samples[sample]['GT'][0] != rec.samples[sample]['GT'][1]:
            yield rec


def get_snp_count


def write_out(out):
    print('ahhh this should be in tsv format', file=out)


def get_args():
    """parse the arguments to this script using argparse"""
    parser = argparse.ArgumentParser(description=
        "Get per-sample counts of reads overlapping STRs "
        "by using counts from nearby SNPs."
    )
    parser.add_argument(
        '-o', '--out', help='output file name (in tsv format)',
        type=argparse.FileType('w'), default=sys.stdout
    )
    parser.add_argument(
        'vcf', help='sorted and phased STRs and SNPs'
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
    variants = get_hets(args.vcf, args.sample)
    # get_variants(args.vcf, args.)
    write_out(args.out)


if __name__ == "__main__":
    main(get_args())
