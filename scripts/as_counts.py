#!/usr/bin/env python
import sys
import argparse



def write_out(fname):
    print('ahhh this should be in tsv format', file=fname)


def get_args():
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
    )
    parser.add_argument(
        'sample', help='the name of the sample for which to extract STR counts'
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    write_out(args.out)


if __name__ == "__main__":
    main()
