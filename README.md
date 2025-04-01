# Snakemake workflow: rna-seq-star-deseq2-GSEA-scInference

This workflow is based on [Snakemake Workflow rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2), but with added functions of GSEA analysis, single cell inference from bulk RNA-seq and time series analysis based on inference data 

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-deseq2).

## Set up before run

Gather necessary information to make modifications of files
config/samples.tsv: sample names, treatment conditions, jointly_handled (to monitor batch effect derived from lab, sequencer etc)
config/units.tsv: sample names, fastq file path, adapter sequence (add -a ADAPTER_SEQUENCE), strandedness (default: none)
config/config.yaml: this is the basic set up configuration. pay special attention to the (diffexp) section. the variables there need to match with config/samples.tsv column header




