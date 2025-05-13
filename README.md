# Snakemake workflow: rna-seq-star-deseq2-GSEA-scInference

This workflow is based on [Snakemake Workflow rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2), but with added functions of GSEA analysis, single cell inference from bulk RNA-seq and time series analysis based on inference data 

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-deseq2)

## Set up before run

Gather necessary information to make modifications of files
config/samples.tsv: sample names, treatment conditions, jointly_handled (to monitor batch effect derived from lab, sequencer etc)

config/units.tsv: sample names, unit names (same sample sequenced in multiple lanes?), fastq file path, adapter sequence (add -a ADAPTER_SEQUENCE), strandedness (default: none)

config/config.yaml: this is the basic set up configuration. the reference files and result directories should be updated for each analysis. pay special attention to the (diffexp) section. the variables there need to match with config/samples.tsv column header

## Workflow structure

The main workflow orchestration file is the main snakemake file called Snakefile. It calls modular snakemake files in in the rules directory, which is a standardized format for the command of each execution step, organized by input, output, shell command or wrapper, logging, parameters, conda run environment pointing to workflow/envs/, which track software versions and should be updated after tool installation.

Custom scripts are saved in scripts directory. They can be called in snakemake file using script directive,  such as "../scripts/deseq2-init.R". This is useful to implement custom python or R scripts, such as deseq analysis.

Instead of pre-downloading tools, many standard tools are hosted in snakemake repository and can be called via wrapper. Only the version of the tools are required in the wrapper directive, such as "v5.10.0/bio/samtools/index". The tools will be automatically downloaded when executing. Current stable repository is v5.10.0 and all wrappers has been updated to the current version.

The sample file data in config.yaml file needs to be processed before being passed to rules at initiation stage, therefore there are input functions in rules/common.smk file to process config/samples.tsv and config/units.tsv files.

config/config.yaml file contains information related to setup of each step of the run. This is different from the rules in that it's used to track sample file data, contrast in diff analysis, whether a specific step is activated or not.

resources directory host downloaded genome fasta, gtf files and indices. When using new genome assembly or annotation files, the following smk files in rules referring to resources need to be modified: align.smk, qc.smk and ref.smk

## Data storage

Snakemake files, associated scripts and resource files are saved in a smaller volume of 50Gb, mounted to EC2 as data1. This volume is backed up with a snapshot image, and can be reused with future projects. Genome index in resources directory needs to be amended if working with a new species.

Config, input and result files are saved in a bigger volume of 500Gb, mounted to EC2 as data2 for processing, and unmounted and deleted after the project is completed and all the data is backed up in s3. Changes need to be made to the files under config directory to capture sample, unit information and contrast customization in deseq analysis.

## Installation

Make sure conda is installed first. You can download it from [here] (https://github.com/conda-forge/miniforge/releases/). For osx, install the x86_64 version because bioconda packages do not support arm64. Then add conda conda/bin to $PATH by modifying ~/.zshrc

```
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/miniforge3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/opt/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/opt/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

conda should be version >=24.7.1. If not

```conda update conda```

Then create the virtual environment

```mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy```

conda-forge and bioconda host necessary packages

Make sure both conda-forge and bioconda channels are available to install necessary packages

```conda config --show channels```

if needed

```conda config --add channels bioconda```
```conda config --add channels conda-forge```

```conda config --set channel_priority strict```

Then navigate to the project directory and activate the environment

```conda activate snakemake```

execute command

```snakemake --cores 2 --sdm conda```

## Issues

1. PackagesNotFoundError: The following packages are not available from current channels:
```- star=2.7.11b*```

Current channels:

```
- https://conda.anaconda.org/conda-forge
- https://conda.anaconda.org/bioconda
```

This is due to my computer running osx-arm64 architecture. Switch to miniforge x86_64 version solved the problem.

linux EC2 instances should not have this issue.

2. STAR alignment error message.

Transcriptome.cpp:18:Transcriptome: exiting because of *INPUT FILE* error: could not open input file /geneInfo.tab

This error has been reported for the latest 2.7.11b version of STAR, and an older version of STAR 2.7.5a worked well according to user feedbacks. So this requires a local installation of the older version (conda deactivate the snakemake environment first)

sudo conda install -c bioconda star=2.7.5a

Then modify align.smk and ref.smk where star is called; change "wrapper" into "shell" so that it will run the locally installed star instead of downloading 2.7.11b from bioconda: /opt/miniforge3/bin/star

Note because I'm running on osx, zcat doesn't work the same way as in linux. The following line need to be changed in align.smk

```--readFilesCommand zcat```

should be changed into

```--readFilesCommand gunzip -c```

## TPM calculation addition
1. Added rule calc_TPM in diffexp.smk

2. Added 
```
- bioconductor-genomicfeatures =1.58
- bioconductor-txdbmaker =1.2
```
to deseq2.yaml dependencies

3. Change get_final_output() in common.smk to reflect TPM count matrix file as part of the final_output

4. Added calc_TPM.R to scripts

## GSEA addition
1. Added rule gsea in diffexp.smk

2. Added gsea.yaml to designate dependencies of bioconda packages and versions compatible with the current operating system:
```
  - bioconductor-org.hs.eg.db =3.20 # this needs to be adjusted based on species used in the experiment
  - r-dplyr =1.1
  - bioconductor-dose =3.6
  - bioconductor-enrichplot =1.26
  - bioconductor-clusterprofiler =4.14
  - r-readr =2.1
  - r-ggridges =0.5
  - r-ggplot2 =3.5
```

3. Change get_final_output() in common.smk to reflect GSEA result file as the final_output

4. Added gsea.R to scripts

5. Pay special attention to library(org.hs.eg.db) in gsea.R. This should correspond to the species used in the experiment

## Single cell inference addition
1. Added rule gsea in sc_infer.smk

2. Added sc_infer.yaml to designate dependencies of bioconda packages and versions compatible with the current operating system:
```
  - bioconductor-granulator =1.14
```

3. Change get_final_output() in common.smk to reflect single cell inference result file as the final_output

4. Added sc_infer.R to scripts

## Need to do

1. test installation and small sample run (complete)

2. add GSEA R script and modify diffexp.smk file to include GSEA rule (complete)

3. add TPM count matrix (complete)

4. add single cell inference R script and smk file (complete)

5. add olink analysis workflow (separate snakemake workflow if no inference needed)

6. full sample run