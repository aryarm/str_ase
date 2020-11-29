# scripts
This directory contains various scripts used by the pipeline.

### [as_counts.py](as_counts.py)
A python script that retrieves unbiased ASE counts of reads overlapping STRs using WASP-corrected counts of reads from SNPs nearby.

### [unphase.py](unphase.py)
A python script that can be used to unphase the genotypes in a VCF in preparation for use by BEAGLE. This script is used by the `prepare` pipeline.

### [remove_empty_alts.py](remove_empty_alts.py)
A python script that can be used to remove variants with empty ALT records in a VCF.
