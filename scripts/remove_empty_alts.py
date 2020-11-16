#!/usr/bin/env python
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])  # auto-detect input format
vcf_out = VariantFile('-', 'w', header=vcf_in.header)

for rec in vcf_in.fetch():
    # don't output variants that have empty alts
    if rec.alts is None:
        continue
    vcf_out.write(rec)
