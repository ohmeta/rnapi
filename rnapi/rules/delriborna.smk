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
        extra = config["params"]["delriborna"]["ribodetector"]["extra"],
        outdir = os.path.join(config["output"]["delriborna"], "short_reads/{sample}")
    log:
        os.path.join(config["output"]["delriborna"], "logs/{sample}.ribodetector.log")
    benchmark:
        os.path.join(config["output"]["delriborna"],
                     "benchmark/ribodetector/{sample}.ribodetector.benchmark.txt")
    threads:
        config["params"]["delriborna"]["threads"]
    conda:
        config["envs"]["delriborna"]
    shell:
        '''
        OUT1={output.reads[0]}
        OUT2={output.reads[1]}
        RNA1={output.reads_rrna[0]}
        RNA2={output.reads_rrna[1]}

        out1=`echo ${{OUT1%.gz}}`
        out2=`echo ${{OUT2%.gz}}`
        rna1=`echo ${{RNA1%.gz}}`
        rna2=`echo ${{RNA2%.gz}}`

        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        {params.ribodetector} \
        --len {params.length} \
        --ensure rrna \
        --input {input.reads} \
        --output $out1 $out2 \
        --rrna $rna1 $rna2 \
        --chunk_size {params.chunk_size} \
        --threads {threads} \
        {params.extra} \
        > {log} 2>&1

        pigz -p {threads} $out1 >> {log} 2>&1
        pigz -p {threads} $out2 >> {log} 2>&1
        pigz -p {threads} $rna1 >> {log} 2>&1
        pigz -p {threads} $rna2 >> {log} 2>&1
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
