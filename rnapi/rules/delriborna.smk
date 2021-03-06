def delriborna_input(wildcards):
    if TRIMMING_DO:
        return get_reads(wildcards, "trimming")
    else:
        return get_reads(wildcards, "raw")


rule delriborna_ribodetector:
    input:
        reads = delriborna_input
    output:
        reads = expand(os.path.join(
            config["output"]["delriborna"],
            "short_reads/{{sample}}/{{sample}}.nonrrna.{read}.fq.gz"),
                       read=["1", "2"]),
        reads_rrna = expand(os.path.join(
            config["output"]["delriborna"],
            "short_reads/{{sample}}/{{sample}}.rrna.{read}.fq.gz"),
                            read=["1", "2"])
    params:
        ribodetector = "ribodetector" \
            if config["params"]["delriborna"]["ribodetector"]["GPU"] \
               else "ribodetector_cpu",
        length = config["params"]["delriborna"]["ribodetector"]["reads_len"],
        chunk_size = config["params"]["delriborna"]["ribodetector"]["chunk_size"],
        extra = config["params"]["delriborna"]["ribodetector"]["extra"]
    benchmark:
        os.path.join(config["output"]["delriborna"],
                     "benchmark/ribodetector/{sample}.ribodetector.tsv")
    threads:
        config["params"]["delriborna"]["threads"]
    log:
        os.path.join(config["output"]["delriborna"],
                     "logs/{sample}.ribodetector.log")
    shell:
        '''
        {params.ribodetector} \
        --len {params.length} \
        --ensure rrna \
        --input {input.reads} \
        --output {output.reads} \
        --rrna {output.reads_rrna} \
        --chunk_size {params.chunk_size} \
        --threads {threads} \
        {params.extra} \
        > {log} 2>&1
        '''


if config["params"]["delriborna"]["ribodetector"]["do"]:
    rule delriborna_ribodetector_all:
        input:
            expand(os.path.join(
                config["output"]["delriborna"],
                "short_reads/{sample}/{sample}.{rna}.{read}.fq.gz"),
                   read=["1", "2"],
                   rna=["nonrrna", "rrna"],
                   sample=SAMPLES.index.unique())

else:
    rule delriborna_ribodetector_all:
        input:


if DELRIBORNA_DO and config["params"]["qcreport"]["do"]:
    rule delriborna_report:
        input:
            lambda wildcards: get_reads(wildcards, "delriborna")
        output:
            os.path.join(config["output"]["delriborna"],
                              "report/stats/{sample}_delriborna_stats.tsv")
        log:
            os.path.join(config["output"]["delriborna"],
                         "logs/report/{sample}.seqkit.log")
        params:
            sample_id = "{sample}",
            fq_encoding = config["params"]["fq_encoding"]
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            shell(
                '''
                seqkit stats \
                --all \
                --basename \
                --tabular \
                --fq-encoding {params.fq_encoding} \
                --out-file {output} \
                --threads {threads} \
                {input} 2> {log}
                ''')

            if IS_PE:
                rnapi.change(output[0], params.sample_id, "delriborna",
                              "pe", ["fq1", "fq2"])
            else:
                rnapi.change(output[0], params.sample_id, "delriborna",
                              "se", ["fq1"])


    rule delriborna_report_merge:
        input:
            expand(
                os.path.join(config["output"]["delriborna"],
                             "report/stats/{sample}_delriborna_stats.tsv"),
                sample=SAMPLES.index.unique())
        output:
            os.path.join(config["output"]["qcreport"], "delriborna_stats.tsv")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            rnapi.merge(input, rnapi.parse, threads, output=output[0])

    rule delriborna_report_all:
        input:
            os.path.join(config["output"]["qcreport"], "delriborna_stats.tsv")

else:
    rule delriborna_report_all:
        input:


rule delriborna_all:
    input:
        rules.delriborna_ribodetector_all.input,
        rules.delriborna_report_all.input
