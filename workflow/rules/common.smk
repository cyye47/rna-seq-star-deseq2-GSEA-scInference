import glob

import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

RESULTS_DIR = config["results_dir"]

def get_final_output():
    final_output = expand(
        "{results_dir}/gsea/{contrast}/gsea_go_dge_all_results.tsv",
        contrast=list(config["diffexp"]["contrasts"].keys()),
        results_dir=RESULTS_DIR
    )
    final_output.append(RESULTS_DIR + "/deseq2/normcounts.symbol.tsv")
    final_output.append(RESULTS_DIR + "/counts/all.symbol.tsv")
    final_output.append(RESULTS_DIR + "/counts/all_tpm.symbol.tsv")
    final_output.append(RESULTS_DIR + "/qc/multiqc_report.html")
    final_output.append(RESULTS_DIR + "/pca.svg")
    final_output.append(RESULTS_DIR + "/sc_infer/cell_proportion.tsv")
    return final_output


validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name, U=unit.unit_name, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name, U=unit.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )


def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
    )
    assert len(files) > 0
    return files


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


def get_fq(wildcards):
    if config["trimming"]["activate"]:
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "{results_dir}/trimmed/{sample}_{unit}_{group}.fastq.gz",
                        group=["R1", "R2"],
                        results_dir=RESULTS_DIR,
                        **wildcards,
                    ),
                )
            )
        # single end sample
        return {
            "fq1": RESULTS_DIR + "/trimmed/{sample}_{unit}_single.fastq.gz".format(**wildcards)
        }
    else:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit)]
        if pd.isna(u["fq1"]):
            # SRA sample (always paired-end for now)
            accession = u["sra"]
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "sra/{accession}_{group}.fastq",
                        accession=accession,
                        group=["R1", "R2"],
                    ),
                )
            )
        if not is_paired_end(wildcards.sample):
            return {"fq1": f"{u.fq1}"}
        else:
            return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}


def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * units.shape[0]


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_bioc_species_name():
    first_letter = config["ref"]["species"][0]
    subspecies = config["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def get_fastqs(wc):
    if config["trimming"]["activate"]:
        return expand(
            "{results_dir}/trimmed/{sample}/{unit}_{read}.fastq.gz",
            unit=units.loc[wc.sample, "unit_name"],
            results_dir=RESULTS_DIR,
            sample=wc.sample,
            read=wc.read,
        )
    unit = units.loc[wc.sample]
    if all(pd.isna(unit["fq1"])):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand(
            "sra/{accession}_{read}.fastq", accession=accession, read=wc.read[-1]
        )
    fq = "fq{}".format(wc.read[-1])
    return units.loc[wc.sample, fq].tolist()


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]
