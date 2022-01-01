#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np


def parse_gene_tab_init(strandedness):
    global STAR_GENE_TAB_COLUMN__
    STAR_GENE_TAB_COLUMN__ = 3

    if pd.isnull(strandedness) or strandedness == "none" or strandedness == "":
        STAR_GENE_TAB_COLUMN__ = 1  # non stranded protocol
        print("Parsing gene tab: using column 2")
    elif strandedness == "forward":
        STAR_GENE_TAB_COLUMN__ = 2  # 3rd column
        print("Parsing gene tab: using column 3")
    elif strandedness == "reverse":
        STAR_GENE_TAB_COLUMN__ = 3  # 4th column, usually for Illumina truseq
        print("Parsing gene tab: using column 4")
    else:
        sys.exit("strandedness is not right")


def parse_gene_tab(genef):
    sample_id = os.path.basename(os.path.dirname(genef))
    if os.path.exists(genef):
        try:
            df = pd.read_csv(
                genef,
                usecols=[0, STAR_GENE_TAB_COLUMN__],
                names=["gene", sample_id],
                skiprows=4,
                sep="\t",
            ).set_index("gene")
        except pd.errors.EmptyDataError:
            print(f"{genef} is empty, please check")
            return None

        if not df.empty:
            return df
        else:
            print(f"{genef} is empty, please check")
            return None
    else:
        sys.exit(f"{genef} is not exists")


def parse_rsem_gene_TPM(genef):
    sample_id = os.path.basename(os.path.dirname(genef))
    if os.path.exists(genef):
        try:
            df = pd.read_csv(
                genef,
                usecols=[0, 5],
                header=0,
                names=["gene_id", sample_id],
                sep="\t",
            ).set_index("gene_id")
        except pd.errors.EmptyDataError:
            print(f"{genef} is empty, please check")
            return None

        if not df.empty:
            return df
        else:
            print(f"{genef} is empty, please check")
            return None
    else:
        sys.exit(f"{genef} is not exists")


def parse_rsem_gene_FPKM(genef):
    sample_id = os.path.basename(os.path.dirname(genef))
    if os.path.exists(genef):
        try:
            df = pd.read_csv(
                genef,
                usecols=[0, 6],
                header=0,
                names=["gene_id", sample_id],
                sep="\t",
            ).set_index("gene_id")
        except pd.errors.EmptyDataError:
            print(f"{genef} is empty, please check")
            return None

        if not df.empty:
            return df
        else:
            print(f"{genef} is empty, please check")
            return None
    else:
        sys.exit(f"{genef} is not exists")


def parse_rsem_transcript_TPM(transcriptf):
    sample_id = os.path.basename(os.path.dirname(transcriptf))
    if os.path.exists(transcriptf):
        try:
            df = pd.read_csv(
                transcriptf,
                usecols=[0, 5],
                header=0,
                names=["transcript_id", sample_id],
                sep="\t",
            ).set_index("transcript_id")
        except pd.errors.EmptyDataError:
            print(f"{transcriptf} is empty, please check")
            return None

        if not df.empty:
            return df
        else:
            print(f"{transcriptf} is empty, please check")
            return None
    else:
        sys.exit(f"{transcriptf} is not exists")


def parse_rsem_transcript_FPKM(transcriptf):
    sample_id = os.path.basename(os.path.dirname(transcriptf))
    if os.path.exists(transcriptf):
        try:
            df = pd.read_csv(
                transcriptf,
                usecols=[0, 6],
                header=0,
                names=["transcript_id", sample_id],
                sep="\t",
            ).set_index("transcript_id")
        except pd.errors.EmptyDataError:
            print(f"{transcriptf} is empty, please check")
            return None

        if not df.empty:
            return df
        else:
            print(f"{transcriptf} is empty, please check")
            return None
    else:
        sys.exit(f"{transcriptf} is not exists")


def parse_salmon_TPM(qf):
    sample_id = os.path.basename(os.path.dirname(qf))
    if os.path.exists(qf):
        try:
            df = pd.read_csv(
                qf,
                usecols=[0, 3],
                header=0,
                names=["id", sample_id],
                sep="\t",
            ).set_index("id")
        except pd.errors.EmptyDataError:
            print(f"{qf} is empty, please check")
            return None

        if not df.empty:
            return df
        else:
            print(f"{qf} is empty, please check")
            return None
    else:
        sys.exit(f"{qf} is not exists")


def parse_salmon_count(qf):
    sample_id = os.path.basename(os.path.dirname(qf))
    if os.path.exists(qf):
        try:
            df = pd.read_csv(
                qf,
                usecols=[0, 4],
                header=0,
                names=["id", sample_id],
                sep="\t",
            ).set_index("id")
        except pd.errors.EmptyDataError:
            print(f"{qf} is empty, please check")
            return None

        if not df.empty:
            return df
        else:
            print(f"{qf} is empty, please check")
            return None
    else:
        sys.exit(f"{qf} is not exists")
