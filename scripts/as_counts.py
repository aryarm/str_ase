#!/usr/bin/env python
import sys
import argparse
import pandas as pd
from pysam import VariantFile
from types import SimpleNamespace



def get_hets(vcf, sample):
    """extract heterozygous STRs and SNPs from the vcf"""
    vcf = VariantFile(vcf)
    vcf.subset_samples([sample])
    # iterate through each variant record in the VCF
    # TODO: make sure you're considering cases where the POS is duplicated
    for rec in vcf.fetch():
        variant = rec.samples[sample]
        is_het = variant['GT'][0] != variant['GT'][1]
        is_snp = not ((len(variant.alleles[0]) - len(variant.alleles[1])))
        # limit our analysis to only variants that are heterozygous
        # and they must be either STRs or SNPs
        if (is_het and (is_snp or rec.id.startswith("STR_"))):
            yield rec


def find_nearest(strs, snps):
    """find the nearest of the SNPs to each STR"""
    for STR in strs:
        # how far is this STR from the last SNP and the current SNP?
        # distance is infinite if they live on different chromosomes
        distances = tuple(
            (
                abs(STR.pos-var[0].pos)
                if var[0].chrom == STR.chrom else float('inf')
            ) for var in snps
        )
        distance = min(distances)
        yield STR, snps[distances.index(distance)], distance


def get_str_snp_pairs(variants, snp_counts):
    """return tuples of (str, nearest snp)"""
    # holds the most recent SNP
    snp = (SimpleNamespace(chrom='', pos=float('inf')), snp_counts[:0])
    distance = float('inf')
    # holds the current STRs
    strs = []
    # iterate through the SNPs and STRs
    for variant in variants:
        # if we found an STR, save it for later
        if variant.id.startswith("STR_"):
            strs.append(variant)
            continue
        # if we're here, it means this variant is a SNP, not an STR
        # first, we must check whether this SNP appears in the snp_counts file
        try:
            snp_count = snp_counts.loc[(variant.chrom, variant.pos)]
        except KeyError:
            continue
        # return pairs of each STR and the nearest SNP
        yield from find_nearest(strs, (snp, (variant, snp_count)))
        strs = []
        # save the current SNP and move on
        snp = variant, snp_count
    # yield any STRs we might have accumulated before reaching the end of the VCF
    if distance != float('inf'):
        yield from find_nearest(strs, (snp,))


def get_snp_counts(snp_counts):
    """parse the SNP ASE table"""
    columns = [
        'CHR', 'POS', 'SAMPLE_ID', 'SUBJECT_ID', 'TISSUE_ID', 'REF_ALLELE',
        'ALT_ALLELE', 'REF_COUNT', 'ALT_COUNT', 'TOTAL_COUNT', 'REF_RATIO',
        'OTHER_ALLELE_COUNT', 'BINOM_P', 'BINOM_P_ADJUSTED', 'GENOTYPE',
        'VARIANT_ANNOTATION', 'GENE_ID'
    ]
    # rename the columns
    # the uppercase is frankly a bit annoying
    new_col_names = {
        col: (col[:-len("_ID")] if col.endswith('_ID') else col).lower()
        for col in columns
    }
    new_col_names['GENOTYPE'] = 'snp_gt'
    new_col_names['REF_ALLELE'] = 'snp_ref'
    new_col_names['ALT_ALLELE'] = 'snp_alt'
    # use pandas to get the datatypes of each column
    types_dict = pd.read_csv(
        snp_counts, sep="\t", header=0, nrows=5, usecols=columns
    ).dtypes
    # set the 'CHR' datatype to np.dtype('O') (aka object)
    types_dict['CHR'] = types_dict['TISSUE_ID']
    # read the entire tsv into memory
    # TODO: read on a stream and convert this function into a generator
    counts = pd.read_csv(
        snp_counts, sep="\t", header=0, dtype=dict(types_dict),
        usecols=columns
    )
    return counts.rename(new_col_names, axis=1).set_index(['chr', 'pos', 'tissue'])


def get_str_counts(pairs, sample):
    """extract the desired information from every STR-SNP pair"""
    for pair in pairs:
        STR, (snp, snp_count), distance = pair
        str_count = snp_count.copy()
        # extract the columns we need:
        # STR ID, gene ID, sample name, the ASE, the location, the STR alleles, the SNP
        # alleles, any other SNP ASE stats, the SNP/STR distance, and the tissue ID
        str_count.insert(0, 'str', STR.id)
        str_count.insert(1, 'snp', snp.id)
        str_count.insert(2, 'chrom', snp.chrom)
        str_count.insert(3, 'str_pos', STR.pos)
        str_count.insert(4, 'snp_pos', snp.pos)
        STR = STR.samples[sample]
        str_count.insert(5, 'str_ref', STR.alleles[0])
        str_count.insert(6, 'str_alt', STR.alleles[1])
        str_count.insert(7, 'distance', distance)
        str_count.insert(8, 'str_gt', "{}|{}".format(*STR['GT']))
        str_count['snp_gt'] = str_count['snp_gt'].str[len('GT;'):]
        yield str_count


def write_out(out, sample, str_counts):
    """write each pd dataframe to the output in append mode"""
    first_line = True
    for str_count in str_counts:
        str_count.to_csv(out, mode='a', sep="\t", header=first_line)
        first_line = False


def get_args():
    """parse the arguments to this script using argparse"""
    parser = argparse.ArgumentParser(description=
        "Get per-sample counts of reads overlapping STRs "
        "by using counts from nearby SNPs."
    )
    parser.add_argument(
        '-o', '--out', help='path to TSV output (default: stdout)',
        type=argparse.FileType('w'), default=sys.stdout
    )
    parser.add_argument(
        'counts', help='WASP-corrected counts of SNPs for this GTEX sample'
        # TODO: describe the structure of this tsv file better
    )
    parser.add_argument(
        'sample', help='the name of the sample for which to extract STR counts'
    )
    parser.add_argument(
        'vcf', default="-", nargs='?',
        help="sorted/phased VCF w/ STRs (prefix ID: 'STR_') and SNPs (default: stdin)"
    )
    args = parser.parse_args()
    return args


def main(args):
    write_out(
        args.out, args.sample,
        get_str_counts(
            get_str_snp_pairs(
                get_hets(args.vcf, args.sample),
                get_snp_counts(args.counts)
            ),
            args.sample
        )
    )


if __name__ == "__main__":
    main(get_args())
