LOG_DIR = config["log_dir"]

rule get_genome:
    output:
        config["ref"]["fa_file"],
    log:
        LOG_DIR + "/get-genome.log",
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
        config["ref"]["gtf_file"],
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        LOG_DIR + "/get_annotation.log",
    wrapper:
        "v5.10.0/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        config["ref"]["fa_file"],
    output:
        config["ref"]["fa.fai_file"],
    log:
        LOG_DIR + "/genome-faidx.log",
    cache: True
    wrapper:
        "v5.10.0/bio/samtools/faidx"


rule bwa_index:
    input:
        config["ref"]["fa_file"],
    output:
        multiext(config["ref"]["fa_file"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        LOG_DIR + "/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v5.10.0/bio/bwa/index"


rule star_index:
    input:
        fasta=config["ref"]["fa_file"],
        gtf=config["ref"]["gtf_file"]
    output:
        directory(config["ref"]["star_dir"])
    threads: 4
    params:
        sjdbOverhang=100,
    log:
        LOG_DIR + "/star_index_genome.log",
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
