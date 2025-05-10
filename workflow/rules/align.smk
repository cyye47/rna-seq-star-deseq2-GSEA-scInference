RESULTS_DIR = config["results_dir"]
LOG_DIR = config["log_dir"]

rule align:
    input:
        unpack(get_fq),
        index=config["ref"]["star_dir"],
        gtf=config["ref"]["gtf_file"]
    output:
        aln = RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene = RESULTS_DIR + "/star/{sample}_{unit}/ReadsPerGene.out.tab"
    threads: 4
    log:
        LOG_DIR + "/star/{sample}_{unit}.log",
    params:
        out_prefix = RESULTS_DIR + "/star/{sample}_{unit}/",
        fq2 = lambda wildcards, input: input.fq2 if "fq2" in input else ""
    shell:
        """
        /opt/miniforge3/bin/star \
        --runThreadN {threads} \
        --genomeDir {input.index} \
        --readFilesIn {input.fq1} {params.fq2}\
        --readFilesCommand gunzip -c \
        --sjdbGTFfile {input.gtf} \
        --outFileNamePrefix {params.out_prefix} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        &> {log}
        """
