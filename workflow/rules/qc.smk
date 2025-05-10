## RSEQC

RESULTS_DIR = config["results_dir"]
LOG_DIR = config["log_dir"]

rule rseqc_gtf2bed:
    input:
        config["ref"]["gtf_file"],
    output:
        bed = RESULTS_DIR + "/qc/rseqc/annotation.bed",
        db=temp(RESULTS_DIR + "/qc/rseqc/annotation.db"),
    log:
        LOG_DIR + "/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam = RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed = RESULTS_DIR + "/qc/rseqc/annotation.bed",
    output:
        RESULTS_DIR + "/qc/rseqc/{sample}_{unit}.junctionanno.junction.bed",
    priority: 1
    log:
        LOG_DIR + "/rseqc/rseqc_junction_annotation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam = RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed = RESULTS_DIR + "/qc/rseqc/annotation.bed",
    output:
        RESULTS_DIR + "/qc/rseqc/{sample}_{unit}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        LOG_DIR + "/rseqc/rseqc_junction_saturation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
    output:
        RESULTS_DIR + "/qc/rseqc/{sample}_{unit}.stats.txt",
    priority: 1
    log:
        LOG_DIR + "/rseqc/rseqc_stat/{sample}_{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam = RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed = RESULTS_DIR + "/qc/rseqc/annotation.bed",
    output:
        RESULTS_DIR + "/qc/rseqc/{sample}_{unit}.infer_experiment.txt",
    priority: 1
    log:
        LOG_DIR + "/rseqc/rseqc_infer/{sample}_{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam = RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed = RESULTS_DIR + "/qc/rseqc/annotation.bed",
    output:
        RESULTS_DIR + "/qc/rseqc/{sample}_{unit}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        LOG_DIR + "/rseqc/rseqc_innerdis/{sample}_{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam = RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed = RESULTS_DIR + "/qc/rseqc/annotation.bed",
    output:
        RESULTS_DIR + "/qc/rseqc/{sample}_{unit}.readdistribution.txt",
    priority: 1
    log:
        LOG_DIR + "/rseqc/rseqc_readdis/{sample}_{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
    output:
        RESULTS_DIR + "/qc/rseqc/{sample}_{unit}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        LOG_DIR + "/rseqc/rseqc_readdup/{sample}_{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        RESULTS_DIR + "/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
    output:
        RESULTS_DIR + "/qc/rseqc/{sample}_{unit}.readgc.GC_plot.pdf",
    priority: 1
    log:
        LOG_DIR + "/rseqc/rseqc_readgc/{sample}_{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule multiqc:
    input:
        expand(
            RESULTS_DIR + "/star/{unit.sample_name}_{unit.unit_name}/Aligned.sortedByCoord.out.bam",
            unit=units.itertuples(),
        ),
        expand(
            RESULTS_DIR + "/qc/rseqc/{unit.sample_name}_{unit.unit_name}.junctionanno.junction.bed",
            unit=units.itertuples(),
        ),
        expand(
            RESULTS_DIR + "/qc/rseqc/{unit.sample_name}_{unit.unit_name}.junctionsat.junctionSaturation_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            RESULTS_DIR + "/qc/rseqc/{unit.sample_name}_{unit.unit_name}.infer_experiment.txt",
            unit=units.itertuples(),
        ),
        expand(
            RESULTS_DIR + "/qc/rseqc/{unit.sample_name}_{unit.unit_name}.stats.txt",
            unit=units.itertuples(),
        ),
        expand(
            RESULTS_DIR + "/qc/rseqc/{unit.sample_name}_{unit.unit_name}.inner_distance_freq.inner_distance.txt",
            unit=units.itertuples(),
        ),
        expand(
            RESULTS_DIR + "/qc/rseqc/{unit.sample_name}_{unit.unit_name}.readdistribution.txt",
            unit=units.itertuples(),
        ),
        expand(
            RESULTS_DIR + "/qc/rseqc/{unit.sample_name}_{unit.unit_name}.readdup.DupRate_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            RESULTS_DIR + "/qc/rseqc/{unit.sample_name}_{unit.unit_name}.readgc.GC_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            LOG_DIR + "/rseqc/rseqc_junction_annotation/{unit.sample_name}_{unit.unit_name}.log",
            unit=units.itertuples(),
        ),
    output:
        RESULTS_DIR + "/qc/multiqc_report.html",
    log:
        LOG_DIR + "/multiqc.log",
    wrapper:
        "v5.10.0/bio/multiqc"
