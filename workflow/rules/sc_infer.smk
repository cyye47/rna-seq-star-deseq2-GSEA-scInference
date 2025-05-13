RESULTS_DIR = config["results_dir"]
LOG_DIR = config["log_dir"]

rule sc_infer:
    input:
        pbmc_reference = config["ref"]["pbmc_reference"],
        tpm = RESULTS_DIR + "/counts/all_tpm.tsv",
    output:
        model_performance = RESULTS_DIR + "/sc_infer/model_performance.pdf",
        proportion_fig = RESULTS_DIR + "/sc_infer/cell_proportion.pdf",
        proportion_table = RESULTS_DIR + "/sc_infer/cell_proportion.tsv",
    log:
        LOG_DIR + "/sc_infer.log"
    conda:
        "../envs/sc_infer.yaml"
    script:
        "../scripts/sc_infer.R"
