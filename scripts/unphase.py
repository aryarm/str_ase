#!/usr/bin/env python
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])  # auto-detect input format
vcf_out = VariantFile('-', 'w', header=vcf_in.header)

for rec in vcf_in.fetch():
    for sample in rec.samples:
        # BEAGLE complains when there's no allele separator
        # so we have to convert '.' to './.'
        if rec.samples[sample]['GT'] == (None,):
            rec.samples[sample]['GT'] = (None, None)
        # set the genotype as unphased
        rec.samples[sample].phased = False
    vcf_out.write(rec)
