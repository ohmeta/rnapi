rule assembly_xcr_trust4:
    input:
        reads = get_clean_reads,
        coordinate = config["params"]["assembly"]["trust4"]["coordinate_fasta"],
        reference = config["params"]["assembly"]["trust4"]["reference_fasta"]
    output:
        expand(os.path.join(config["output"]["assembly"],
                            "trust4/{{sample}}/{{sample}}_{suffix}"),
               suffix=["airr.tsv", "annot.fa", "assembled_reads.fa",
                       "cdr3.out", "final.out", "raw.out", "report.tsv",
                       "toassemble_1.fq", "toassemble_2.fq"])
    threads:
        config["params"]["assembly"]["threads"]
    benchmark:
        os.path.join(config["output"]["assembly"],
                     "benchmark/trust4/{sample}.trust4.benchmark.txt")
    log:
        os.path.join(config["output"]["assembly"],
                     "logs/trust4/{sample}.trust4.log")
    params:
        outdir = os.path.join(config["output"]["assembly"], "trust4/{sample}"),
        prefix = "{sample}"
    shell:
        '''
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        run-trust4 \
        -f {input.coordinate} \
        --ref {input.reference} \
        -1 {input.reads[0]} \
        -2 {input.reads[1]} \
        -o {params.prefix} \
        --od {params.outdir} \
        -t {threads} \
        --stage 0 \
        >{log} 2>&1
        '''


if config["params"]["assembly"]["trust4"]["do"]:
    rule assembly_xcr_trust4_all:
        input:
            expand(os.path.join(config["output"]["assembly"],
                            "trust4/{sample}/{sample}_{suffix}"),
                   suffix=["airr.tsv", "annot.fa", "assembled_reads.fa",
                           "cdr3.out", "final.out", "raw.out", "report.tsv",
                           "toassemble_1.fq", "toassemble_2.fq"],
                   sample=SAMPLES.index.unique())

else:
    rule assembly_xcr_trust4_all:
        input:


rule assembly_all:
    input:
        rules.assembly_xcr_trust4_all.input
