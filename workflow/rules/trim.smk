RESULTS_DIR = config["results_dir"]
LOG_DIR = config["log_dir"]

rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    log:
        LOG_DIR + "/get-sra/{accession}.log",
    wrapper:
        "v5.10.0/bio/sra-tools/fasterq-dump"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{unit}.{fq}.{ext}"),
    log:
        LOG_DIR + "/pipe-fastqs/catadapt/{sample}_{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1 = RESULTS_DIR + "/trimmed/{sample}_{unit}_R1.fastq.gz",
        fastq2 = RESULTS_DIR + "/trimmed/{sample}_{unit}_R2.fastq.gz",
        qc = RESULTS_DIR + "/trimmed/{sample}_{unit}.paired.qc.txt",
    log:
        LOG_DIR + "/cutadapt/{sample}_{unit}.log",
    params:
        extra=config["params"]["cutadapt-pe"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "v5.10.0/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq = RESULTS_DIR + "/trimmed/{sample}_{unit}_single.fastq.gz",
        qc = RESULTS_DIR + "/trimmed/{sample}_{unit}_single.qc.txt",
    log:
        LOG_DIR + "/cutadapt/{sample}_{unit}.log",
    params:
        extra=config["params"]["cutadapt-se"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "v5.10.0/bio/cutadapt/se"
