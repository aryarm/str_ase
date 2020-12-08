# scripts
This directory contains various scripts used by the pipeline.

### [as_counts.py](as_counts.py)
A python script that retrieves unbiased ASE counts of reads overlapping STRs using WASP-corrected counts of reads from SNPs nearby.

### [liftover_counts.py](liftover_counts.py)
A python script for lifting over TSVs whose first two columns are CHROM and POS, respectively.

### [plot_associations.py](plot_associations.py)
A python script to create plots of allele specific expression vs STR repeat number. It also fits simple linear regression models to the plots.

### [prioritize_STRs.py](prioritize_STRs.py)
A python script for sorting STRs by a custom importance metric consisting of the number of supporting samples, the distance from the nearest SNP, and the strength of the ASE at the nearest SNP.

### [remove_empty_alts.py](remove_empty_alts.py)
A python script that can be used to remove variants with empty ALT records in a VCF.

### [unphase.py](unphase.py)
A python script that can be used to unphase the genotypes in a VCF in preparation for use by BEAGLE. This script is used by the `prepare` pipeline.
