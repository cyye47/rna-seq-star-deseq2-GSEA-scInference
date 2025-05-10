RESULTS_DIR = config["results_dir"]
LOG_DIR = config["log_dir"]

rule count_matrix:
    input:
        expand(
            RESULTS_DIR + "/star/{unit.sample_name}_{unit.unit_name}/ReadsPerGene.out.tab",
            unit=units.itertuples()
        ),
    output:
        RESULTS_DIR + "/counts/all.tsv",
    log:
        LOG_DIR + "/count-matrix.log",
    params:
        samples=units["sample_name"].tolist(),
        strand=get_strandedness(units),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule gene_2_symbol:
    input:
        counts="{prefix}.tsv",
    output:
        symbol="{prefix}.symbol.tsv",
    params:
        species=get_bioc_species_name(),
    #log:
        #LOG_DIR + "/gene2symbol/{prefix}.log",
    conda:
        "../envs/biomart.yaml"
    script:
        "../scripts/gene2symbol.R"


rule deseq2_init:
    input:
        counts = RESULTS_DIR + "/counts/all.tsv",
    output:
        RESULTS_DIR + "/deseq2/all.rds",
        RESULTS_DIR + "/deseq2/normcounts.tsv",
    conda:
        "../envs/deseq2.yaml"
    log:
        LOG_DIR + "/deseq2/init.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        RESULTS_DIR + "/deseq2/all.rds",
    output:
        report(RESULTS_DIR + "/pca.{variable}.svg", "../report/pca.rst"),
    conda:
        "../envs/deseq2.yaml"
    log:
        LOG_DIR + "/pca.{variable}.log",
    script:
        "../scripts/plot-pca.R"


rule deseq2:
    input:
        RESULTS_DIR + "/deseq2/all.rds",
    output:
        table=report(RESULTS_DIR + "/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report(RESULTS_DIR + "/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        LOG_DIR + "/deseq2/{contrast}.diffexp.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2.R"

rule gsea:
    input:
        diffexpr = RESULTS_DIR + "/diffexp/{contrast}.diffexp.tsv"
    output:
        output_go = RESULTS_DIR + "/gsea/{contrast}/gsea_go_dge_all_results.tsv",
        output_dir = directory(RESULTS_DIR + "/gsea/{contrast}")
    params:
        contrast=get_contrast,
    conda:
        "../envs/gsea.yaml"
    log:
        LOG_DIR + "/gsea/{contrast}.gsea.log",
    script:
        "../scripts/gsea.R"