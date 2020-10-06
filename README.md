[![Snakemake](https://img.shields.io/badge/snakemake-≥5.26.1-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# str_ase
Allele-specific expression (ASE) for short tandem repeats (STRs)

# download
Execute the following command.
```
git clone https://github.com/aryarm/str_ase
```

# setup
The pipeline is written as a Snakefile which can be executed via [Snakemake](https://snakemake.readthedocs.io). We recommend installing version 5.26.1:
```
conda create -n snakemake -c conda-forge --no-channel-priority 'bioconda::snakemake==5.26.1'
```
We highly recommend you install [Snakemake via conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda) like this so that you can use the `--use-conda` flag when calling `snakemake` to let it [automatically handle all dependencies](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) of the pipeline. Otherwise, you must manually install the dependencies listed in the [env file](env.yml).

# execution
1. Activate snakemake via `conda`:
    ```
    conda activate snakemake
    ```
2. Execute the pipeline

    Locally:
    ```
    ./run &
    ```
Log files describing the output of the pipeline will be created within the output directory. The `log` file contains a basic description of the progress of each rule, while the `qlog` file is more detailed.

### Executing the pipeline on your own data
You must modify [the config.yaml file](config.yml) to specify paths to your data before you perform step 2 above.

### If this is your first time using Snakemake
We recommend that you run `snakemake --help` to learn about Snakemake's options. For example, to check that the pipeline will be executed correctly before you run it, you can call Snakemake with the `-n -p -r` flags. This is also a good way to familiarize yourself with the steps of the pipeline and their inputs and outputs (the latter of which are inputs to the first rule in each workflow -- ie the `all` rule).

Note that Snakemake will not recreate output that it has already generated, unless you request it. If a job fails or is interrupted, subsequent executions of Snakemake will just pick up where it left off. This can also apply to files that *you* create and provide in place of the files it would have generated.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a BAM). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipeline.

# files
### Snakefile
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline defining rules for every step of the analysis.

### config.yml
Defines options and input for the Snakemake pipeline.

### run
An example bash script for executing the pipeline using `snakemake` and `conda`. Any arguments to this script are passed directly to `snakemake`.

