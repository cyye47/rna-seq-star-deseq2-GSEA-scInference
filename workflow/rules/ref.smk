rule get_genome:
    output:
        "/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genome.fa",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "v5.10.0/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genes.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "logs/get_annotation.log",
    wrapper:
        "v5.10.0/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        "/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genome.fa",
    output:
        "/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genome.fa.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v5.10.0/bio/samtools/faidx"


rule bwa_index:
    input:
        "/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genome.fa",
    output:
        multiext("/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v5.10.0/bio/bwa/index"


rule star_index:
    input:
        fasta="/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genome.fa",
        gtf="/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/genes.gtf"
    output:
        directory("/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/scerevisiae_R64-1-1/star_genome")
    threads: 4
    params:
        sjdbOverhang=100,
    log:
        "logs/star_index_genome.log",
    cache: True
    shell:
        """
        /opt/miniforge3/bin/star \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.sjdbOverhang} \
        &> {log}
        """
