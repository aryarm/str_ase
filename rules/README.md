# Snakemake rules

### prepare.smk
Prepares a small, example dataset for use by the pipeline. This is stored in the `data/` directory.
#### configuration
The `prepare` pipeline uses the `rules/prepare.yml` and `config.yml` files for configuration. You must fill out these config files prior to executing the `prepare` pipeline.
#### execution
You can execute this file on its own, as its own pipeline, or extend the main pipeline with the rules from this one. To execute it on its own, do
```
./run -s rules/prepare.smk &
```
