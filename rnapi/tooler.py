#!/usr/bin/env python3

import os
import concurrent.futures
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json

import matplotlib

matplotlib.use("agg")


def parse(stats_file):
    if os.path.exists(stats_file):
        try:
            df = pd.read_csv(stats_file, sep="\t")
        except pd.errors.EmptyDataError:
            print("%s is empty, please check" % stats_file)
            return None

        if not df.empty:
            return df
        else:
            return None
    else:
        print("%s is not exists" % stats_file)
        return None


def merge(input_list, func, workers, **kwargs):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df in executor.map(func, input_list):
            if df is not None:
                df_list.append(df)

    df_ = pd.concat(df_list)

    if "output" in kwargs:
        df_.to_csv(kwargs["output"], sep="\t", index=False)
    return df_


def merge_cols(input_list, func, workers, **kwargs):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df in executor.map(func, input_list):
            if df is not None:
                df_list.append(df)

    df_ = pd.concat(df_list, axis=1)

    if "output" in kwargs:
        df_.to_csv(kwargs["output"], sep="\t", index=False)
    return df_


def merge2(input_list, func, workers, **kwargs):
    df1_list = []
    df2_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df1, df2 in executor.map(func, input_list):
            if df1 is not None:
                df1_list.append(df1)
            if df2 is not None:
                df2_list.append(df2)

    df_1 = pd.concat(df1_list)
    df_2 = pd.concat(df2_list)

    if "output_1" in kwargs:
        df_1.to_csv(kwargs["output_1"], sep="\t", index=False)
    if "output_2" in kwargs:
        df_2.to_csv(kwargs["output_2"], sep="\t", index=False)

    return df_1, df_2


def parse_fastp_json(json_f, paired):
    trimming_dict = {}
    sample_id = os.path.basename(json_f).split(".")[0]
    with open(json_f, "r") as ih:
        json_ob = json.load(ih)
        trimming_dict["sample_id"] = sample_id
        trimming_dict["before_filtering_total_reads"] = json_ob["summary"][
            "before_filtering"
        ]["total_reads"]
        trimming_dict["before_filtering_total_bases"] = json_ob["summary"][
            "before_filtering"
        ]["total_bases"]
        trimming_dict["before_filtering_q20_bases"] = json_ob["summary"][
            "before_filtering"
        ]["q20_bases"]
        trimming_dict["before_filtering_q30_bases"] = json_ob["summary"][
            "before_filtering"
        ]["q30_bases"]
        trimming_dict["before_filtering_q20_rate"] = json_ob["summary"][
            "before_filtering"
        ]["q20_rate"]
        trimming_dict["before_filtering_q30_rate"] = json_ob["summary"][
            "before_filtering"
        ]["q30_rate"]
        trimming_dict["before_filtering_qc_content"] = json_ob["summary"][
            "before_filtering"
        ]["qc_content"]
        trimming_dict["before_filtering_read1_mean_length"] = json_ob["summary"][
            "before_filtering"
        ]["read1_mean_length"]
        if paired:
            trimming_dict["before_filtering_read2_mean_length"] = json_ob["summary"][
                "before_filtering"
            ]["read2_mean_length"]

        trimming_dict["read1_before_filtering_total_reads"] = json_ob[
            "read1_before_filtering"
        ]["total_reads"]
        trimming_dict["read1_before_filtering_total_bases"] = json_ob[
            "read1_before_filtering"
        ]["total_bases"]
        trimming_dict["read1_before_filtering_q20_bases"] = json_ob[
            "read1_before_filtering"
        ]["q20_bases"]
        trimming_dict["read1_before_filtering_q30_bases"] = json_ob[
            "read1_before_filtering"
        ]["q30_bases"]
        if paired:
            trimming_dict["read2_before_filtering_total_reads"] = json_ob[
                "read2_before_filtering"
            ]["total_reads"]
            trimming_dict["read2_before_filtering_total_bases"] = json_ob[
                "read2_before_filtering"
            ]["total_bases"]
            trimming_dict["read2_before_filtering_q20_bases"] = json_ob[
                "read2_before_filtering"
            ]["q20_bases"]
            trimming_dict["read2_before_filtering_q30_bases"] = json_ob[
                "read2_before_filtering"
            ]["q30_bases"]

        trimming_dict["after_filtering_total_reads"] = json_ob["summary"][
            "after_filtering"
        ]["total_reads"]
        trimming_dict["after_filtering_total_bases"] = json_ob["summary"][
            "after_filtering"
        ]["total_bases"]
        trimming_dict["after_filtering_q20_bases"] = json_ob["summary"][
            "after_filtering"
        ]["q20_bases"]
        trimming_dict["after_filtering_q30_bases"] = json_ob["summary"][
            "after_filtering"
        ]["q30_bases"]
        trimming_dict["after_filtering_q20_rate"] = json_ob["summary"][
            "after_filtering"
        ]["q20_rate"]
        trimming_dict["after_filtering_q30_rate"] = json_ob["summary"][
            "after_filtering"
        ]["q30_rate"]
        trimming_dict["after_filtering_qc_content"] = json_ob["summary"][
            "after_filtering"
        ]["qc_content"]
        trimming_dict["after_filtering_read1_mean_length"] = json_ob["summary"][
            "after_filtering"
        ]["read1_mean_length"]
        if paired:
            trimming_dict["after_filtering_read2_mean_length"] = json_ob["summary"][
                "after_filtering"
            ]["read2_mean_length"]

        trimming_dict["read1_after_filtering_total_reads"] = json_ob[
            "read1_after_filtering"
        ]["total_reads"]
        trimming_dict["read1_after_filtering_total_bases"] = json_ob[
            "read1_after_filtering"
        ]["total_bases"]
        trimming_dict["read1_after_filtering_q20_bases"] = json_ob[
            "read1_after_filtering"
        ]["q20_bases"]
        trimming_dict["read1_after_filtering_q30_bases"] = json_ob[
            "read1_after_filtering"
        ]["q30_bases"]
        if paired:
            trimming_dict["read2_after_filtering_total_reads"] = json_ob[
                "read2_after_filtering"
            ]["total_reads"]
            trimming_dict["read2_after_filtering_total_bases"] = json_ob[
                "read2_after_filtering"
            ]["total_bases"]
            trimming_dict["read2_after_filtering_q20_bases"] = json_ob[
                "read2_after_filtering"
            ]["q20_bases"]
            trimming_dict["read2_after_filtering_q30_bases"] = json_ob[
                "read2_after_filtering"
            ]["q30_bases"]

        return trimming_dict


def change(output, sample_id, step, fq_type, reads_list):
    df = pd.read_csv(output, sep="\t").sort_values("file", ascending=True)
    df["id"] = sample_id
    df["reads"] = reads_list
    df["step"] = step
    df["fq_type"] = fq_type
    df.to_csv(output, sep="\t", index=False)


def qc_bar_plot(df, engine, stacked=False, **kwargs):
    if engine == "seaborn":
        # seaborn don't like stacked barplot
        f, ax = plt.subplots(figsize=(10, 7))
        df_ = df.query('reads=="fq1"')
        sns.barplot(x="id", y="num_seqs", hue="step", data=df_)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=-90)

    elif engine == "pandas":
        if not stacked:
            df_ = (
                df.query('reads=="fq1"')
                .pivot(index="id", columns="step", values="num_seqs")
                .loc[:, ["raw", "trimming", "rmhost"]]
            )
            df_.plot(kind="bar", figsize=(10, 7))

        else:
            dict_ = {"id": [], "clean": [], "rmhost": [], "trim": []}
            df = df.set_index("id")
            for i in df.index.unique():
                reads_total = 0

                reads_trimmed = 0
                reads_host = 0
                reads_clean = 0

                if (
                    not df.loc[
                        [i],
                    ]
                    .query('reads=="fq1" and step=="raw"')
                    .empty
                ):
                    reads_total = df.loc[[i],].query('reads=="fq1" and step=="raw"')[
                        "num_seqs"
                    ][0]

                if (
                    not df.loc[
                        [i],
                    ]
                    .query('reads=="fq1" and step=="trimming"')
                    .empty
                ):
                    reads_trim = df.loc[[i],].query(
                        'reads=="fq1" and step=="trimming"'
                    )["num_seqs"][0]

                if (
                    not df.loc[
                        [i],
                    ]
                    .query('reads=="fq1" and step=="rmhost"')
                    .empty
                ):
                    reads_clean = df.loc[[i],].query('reads=="fq1" and step=="rmhost"')[
                        "num_seqs"
                    ][0]

                reads_trimmed = reads_total - reads_trim
                reads_host = reads_trim - reads_clean

                dict_["id"].append(i)
                dict_["trim"].append(reads_trimmed)
                dict_["rmhost"].append(reads_host)
                dict_["clean"].append(reads_clean)

            df_ = pd.DataFrame(dict_).sort_values("id").set_index("id")

            colors = ["#2ca02c", "#ff7f0e", "#1f77b4"]

            df_.plot(kind="bar", stacked=True, color=colors, figsize=(10, 7))

    plt.xlabel("Sample ID")
    plt.ylabel("The number of reads(-pair)")
    plt.title("Fastq quality control barplot", fontsize=11)

    if "output" in kwargs:
        plt.savefig(kwargs["output"])
