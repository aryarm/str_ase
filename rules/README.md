# Snakemake rules

### prepare.smk
Prepares a dataset for use by the pipeline. This file can be executed on its own, as its own pipeline. It uses `prepare.yml` for configuration. It can be executed like this
```
./run -s rules/prepare.smk &
```
