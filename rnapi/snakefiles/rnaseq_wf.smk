#!/usr/bin/env snakemake

import sys
import rnapi
import pandas as pd

shell.executable("bash")


IS_PE = True \
    if config["params"]["reads_layout"] == "pe" \
       else False


IS_INTERLEAVED = True \
    if config["params"]["interleaved"] \
       else False


TRIMMING_DO = True \
    if config["params"]["trimming"]["fastp"]["do"] \
       else False


SAMPLES = rnapi.parse_samples(config["params"]["samples"],
                              config["params"]["interleaved"],
                              config["params"]["reads_layout"],
                              config["params"]["begin"])


READS_FORMAT = "sra" \
    if "sra" in SAMPLES.columns \
       else "fastq"


include: "../rules/raw.smk"
include: "../rules/trimming.smk"
include: "../rules/index.smk"
include: "../rules/align.smk"
include: "../rules/quantify.smk"


rule all:
    input:
        rules.raw_all.input,
        rules.trimming_all.input,
        rules.align_all.input,
        rules.quantify_all.input


localrules: \
    prepare_short_reads_all, \
    prepare_reads_all, \
    raw_fastqc_all, \
    raw_report_all, \
    trimming_fastp_all, \
    trimming_report_all, \
    trimming_all, \
    align_genome_star_all, \
    align_transcriptome_star_all, \
    align_star_all, \
    align_all, \
    quantify_genome_star_all, \
    quantify_transcriptome_star_all, \
    quantify_all
