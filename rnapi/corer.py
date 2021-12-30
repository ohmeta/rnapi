#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import textwrap
from io import StringIO

import pandas as pd
import numpy as np

import rnapi

RNASEQ_WF = [
    "prepare_short_reads_all",
    "raw_fastqc_all",
    "raw_report_all",
    "raw_all",
    "trimming_fastp_all",
    "trimming_report_all",
    "trimming_all",
    "align_genome_star_all",
    "align_transcriptome_star_all",
    "align_star_all",
    "align_hisat2_all",
    "align_all",
    "quantify_genome_star_all",
    "quantify_all",
    "pseudo_align_salmon_all",
    "pseudo_align_kallisto_all",
    "quantification_salmon_all",
    "quantification_sleuth_all",
    "all"
]

SCRNASEQ_WF = ["all"]


def run_snakemake(args, unknown, snakefile, workflow):
    conf = rnapi.parse_yaml(args.config)

    if not os.path.exists(conf["params"]["samples"]):
        print("Please specific samples list on init step or change config.yaml manualy")
        sys.exit(1)

    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        args.config
    ] + unknown

    if args.conda_create_envs_only:
        cmd += ["--use-conda", "--conda-create-envs-only"]
    else:
        cmd += [
            "--rerun-incomplete",
            "--keep-going",
            "--printshellcmds",
            "--reason",
            "--until",
            args.task
        ]

        if args.use_conda:
            cmd += ["--use-conda"]
            if args.conda_prefix is not None:
                cmd += ["--conda-prefix", args.conda_prefix]

        if args.list:
            cmd += ["--list"]
        elif args.run_local:
            cmd += ["--cores", str(args.cores)]
        elif args.run_remote:
            cmd += ["--profile", args.profile, "--local-cores", str(args.local_cores), "--jobs", str(args.jobs)]
        elif args.debug:
            cmd += ["--debug-dag", "--dry-run"]
        elif args.dry_run:
            cmd += ["--dry-run"]

    cmd_str = " ".join(cmd).strip()
    print("Running rnapi %s:\n%s" % (workflow, cmd_str))

    env = os.environ.copy()
    proc = subprocess.Popen(
        cmd_str,
        shell=True,
        stdout=sys.stdout,
        stderr=sys.stderr,
        env=env,
    )
    proc.communicate()


def init(args, unknown):
    if args.workdir:
        project = rnapi.rnaconfig(args.workdir)
        print(project.__str__())
        project.create_dirs()
        conf = project.get_config()

        for env_name in conf["envs"]:
            conf["envs"][env_name] = os.path.join(
                os.path.realpath(args.workdir), f"envs/{env_name}.yaml"
            )

        if args.samples:
            conf["params"]["samples"] = args.samples
        else:
            print("Please supply samples table")
            sys.exit(-1)

        rnapi.update_config(
            project.config_file, project.new_config_file, conf, remove=False
        )
    else:
        print("Please supply a workdir!")
        sys.exit(-1)


def rnaseq_wf(args, unknown):
    snakefile = os.path.join(os.path.dirname(__file__), "snakefiles/rnaseq_wf.smk")
    run_snakemake(args, unknown, snakefile, "rnaseq_wf")


def scrnaseq_wf(args, unknown):
    snakefile = os.path.join(os.path.dirname(__file__), "snakefiles/scrnaseq_wf.smk")
    run_snakemake(args, unknown, snakefile, "rnaseq_wf")


def snakemake_summary(snakefile, configfile, task):
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        configfile,
        "--until",
        task,
        "--summary",
    ]
    cmd_out = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    summary = pd.read_csv(StringIO(cmd_out.stdout.read().decode()), sep="\t")
    return summary


def main():
    banner = """

    ██████╗░███╗░░██╗░█████╗░██████╗░██╗
    ██╔══██╗████╗░██║██╔══██╗██╔══██╗██║
    ██████╔╝██╔██╗██║███████║██████╔╝██║
    ██╔══██╗██║╚████║██╔══██║██╔═══╝░██║
    ██║░░██║██║░╚███║██║░░██║██║░░░░░██║
    ╚═╝░░╚═╝╚═╝░░╚══╝╚═╝░░╚═╝╚═╝░░░░░╚═╝

    Omics for All, Open Source for All

      RNA sequence analysis pipeline

"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(banner),
        prog="rnapi",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit",
    )

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "-d",
        "--workdir",
        metavar="WORKDIR",
        type=str,
        default="./",
        help="project workdir",
    )
    common_parser.add_argument(
        "--check-samples",
        dest="check_samples",
        default=False,
        action="store_true",
        help="check samples, default: False",
    )

    run_parser = argparse.ArgumentParser(add_help=False)
    run_parser.add_argument(
        "--config",
        type=str,
        default="./config.yaml",
        help="config.yaml",
    )
    run_parser.add_argument(
        "--profile",
        type=str,
        default="./profiles/slurm",
        help="cluster profile name",
    )
    run_parser.add_argument(
        "--cores",
        type=int,
        default=32,
        help="all job cores, available on '--run-local'")
    run_parser.add_argument(
        "--local-cores",
        type=int,
        dest="local_cores",
        default=8,
        help="local job cores, available on '--run-remote'")
    run_parser.add_argument(
        "--jobs",
        type=int,
        default=80,
        help="cluster job numbers, available on '--run-remote'")
    run_parser.add_argument(
        "--list",
        default=False,
        action="store_true",
        help="list pipeline rules",
    )
    run_parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="debug pipeline",
    )
    run_parser.add_argument(
        "--dry-run",
        default=True,
        dest="dry_run",
        action="store_true",
        help="dry run pipeline",
    )
    run_parser.add_argument(
        "--run-local",
        default=False,
        action="store_true",
        help="run pipeline on local computer",
    )
    run_parser.add_argument(
        "--run-remote",
        default=False,
        action="store_true",
        help="run pipeline on remote cluster",
    )
    run_parser.add_argument(
        "--cluster-engine",
        default="slurm",
        choices=["slurm", "sge", "lsf", "pbs-torque"],
        help="cluster workflow manager engine, support slurm(sbatch) and sge(qsub)"
    )
    run_parser.add_argument("--wait", type=int, default=60, help="wait given seconds")
    run_parser.add_argument(
        "--use-conda",
        default=False,
        dest="use_conda",
        action="store_true",
        help="use conda environment",
    )
    run_parser.add_argument(
        "--conda-prefix",
        default="~/.conda/envs",
        dest="conda_prefix",
        help="conda environment prefix",
    )
    run_parser.add_argument(
        "--conda-create-envs-only",
        default=False,
        dest="conda_create_envs_only",
        action="store_true",
        help="conda create environments only",
    )

    subparsers = parser.add_subparsers(title="available subcommands", metavar="")
    parser_init = subparsers.add_parser(
        "init",
        formatter_class=rnapi.custom_help_formatter,
        parents=[common_parser],
        prog="rnapi init",
        help="init project",
    )
    parser_rnaseq_wf = subparsers.add_parser(
        "rnaseq_wf",
        formatter_class=rnapi.custom_help_formatter,
        parents=[common_parser, run_parser],
        prog="rnapi rnaseq_wf",
        help="RNA seq analysis pipeline",
    )
    parser_scrnaseq_wf = subparsers.add_parser(
        "scrnaseq_wf",
        formatter_class=rnapi.custom_help_formatter,
        parents=[common_parser, run_parser],
        prog="rnapi scrnaseq_wf",
        help="scRNA seq analysis pipeline",
    )

    parser_init.add_argument(
        "-s",
        "--samples",
        type=str,
        default=None,
        help="""desired input:
samples list, tsv format required.
    if it is fastq:
        the header is: [id, fq1, fq2]
    if it is sra:
        the header is: [id, sra]
""",
    )
    parser_init.set_defaults(func=init)

    parser_rnaseq_wf.add_argument(
        "task",
        metavar="TASK",
        nargs="?",
        type=str,
        default="all",
        choices=RNASEQ_WF,
        help="pipeline end point. Allowed values are " + ", ".join(RNASEQ_WF),
    )
    parser_rnaseq_wf.set_defaults(func=rnaseq_wf)

    parser_scrnaseq_wf.add_argument(
        "task",
        metavar="TASK",
        nargs="?",
        type=str,
        default="all",
        choices=SCRNASEQ_WF,
        help="pipeline end point. Allowed values are " + ", ".join(SCRNASEQ_WF),
    )
    parser_scrnaseq_wf.set_defaults(func=scrnaseq_wf)


    args, unknown = parser.parse_known_args()

    try:
        if args.version:
            print("rnapi version %s" % rnapi.__version__)
            sys.exit(0)
        args.func(args, unknown)
    except AttributeError as e:
        print(e)
        parser.print_help()


if __name__ == "__main__":
    main()
