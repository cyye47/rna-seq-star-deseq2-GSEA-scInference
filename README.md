# Snakemake workflow: rna-seq-star-deseq2-GSEA-scInference

This workflow is based on [Snakemake Workflow rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2), but with added functions of GSEA analysis, single cell inference from bulk RNA-seq and time series analysis based on inference data 

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-deseq2).

## Set up before run

Gather necessary information to make modifications of files
config/samples.tsv: sample names, treatment conditions, jointly_handled (to monitor batch effect derived from lab, sequencer etc)

config/units.tsv: sample names, unit names (same sample sequenced in multiple lanes?), fastq file path, adapter sequence (add -a ADAPTER_SEQUENCE), strandedness (default: none)

config/config.yaml: this is the basic set up configuration. pay special attention to the (diffexp) section. the variables there need to match with config/samples.tsv column header

## Workflow structure

The main workflow orchestration file is the main snakemake file called Snakefile. It calls modular snakemake files in in the rules directory, which is a standardized format for the command of each execution step, organized by input, output, shell command or wrapper, logging, parameters, conda run environment pointing to workflow/envs/, which will track software versions and should be updated after tool installation

Custom scripts are saved in scripts directory. They can be called in snakemake file using script directive,  such as "../scripts/deseq2-init.R". This is useful to implement custom python or R scripts, such as deseq analysis.

Instead of pre-downloading tools, many standard tools are hosted in snakemake repository and can be called via wrapper. Only the version of the tools are required in the wrapper directive, such as "v5.10.0/bio/samtools/index". The tools will be automatically downloaded when executing. Current stable repository is v5.10.0 and all wrappers has been updated to the current version.

The sample file data in config.yaml file needs to be processed before being passed to rules at initiation stage, therefore there are input functions in rules/common.smk file to process config/samples.tsv and config/units.tsv files

config/config.yaml file contains information related to setup of each step of the run. This is different from the rules in that it's used to track sample file data, contrast in diff analysis, whether a specific step is activated or not