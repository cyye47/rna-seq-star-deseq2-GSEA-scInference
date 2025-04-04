rule align:
    input:
        unpack(get_fq),
        index="/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/star_genome",
        gtf="/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genes.gtf"
    output:
        aln="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}_{unit}/ReadsPerGene.out.tab"
    threads: 4
    log:
        "logs/star/{sample}_{unit}.log",
    params:
        out_prefix="results/star/{sample}_{unit}/",
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
