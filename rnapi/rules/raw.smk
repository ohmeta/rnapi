def raw_short_reads(wildcards):
    if READS_FORMAT == "fastq":
        if IS_PE:
            if IS_INTERLEAVED:
                return [rnapi.get_reads(SAMPLES, wildcards, "fq1")]
            else:
                return [rnapi.get_reads(SAMPLES, wildcards, "fq1"),
                        rnapi.get_reads(SAMPLES, wildcards, "fq2")]
        else:
            return [rnapi.get_reads(SAMPLES, wildcards, "fq1")]
    elif READS_FORMAT == "sra":
        return [rnapi.get_reads(SAMPLES, wildcards, "sra")]


def short_reads_suffix():
    if IS_PE:
        return [".1", ".2"]
    else:
        return [""]


rule prepare_short_reads:
    input:
        unpack(raw_short_reads)
    output:
        reads = expand(
            os.path.join(
                config["output"]["raw"],
                "short_reads/{{sample}}/{{sample}}.raw{read}.fq.gz"),
            read=short_reads_suffix()) \
            if config["params"]["raw"]["save_reads"] else \
            temp(expand(
                os.path.join(
                    config["output"]["raw"],
                    "short_reads/{{sample}}/{{sample}}.raw{read}.fq.gz"),
                read=short_reads_suffix()))
    params:
        output_dir = os.path.join(config["output"]["raw"], "short_reads/{sample}"),
        interleaved = config["params"]["interleaved"]
    threads:
        config["params"]["raw"]["threads"]
    log:
        os.path.join(config["output"]["raw"], "logs/{sample}_prepare.log")
    run:
        reads_num = len(input)

        if READS_FORMAT == "fastq":
            if IS_PE:
                if not params.interleaved:
                    if reads_num == 2:
                        os.symlink(os.path.realpath(input[0]), output.reads[0])
                        os.symlink(os.path.realpath(input[1]), output.reads[1])
                    else:
                        shell('''cat %s > %s''' % (" ".join(input[0:reads_num//2]), output.reads[0]))
                        shell('''cat %s > %s''' % (" ".join(input[reads_num//2:]), output.reads[1]))
                else:
                    shell(
                        '''
                        cat {input} | \
                        tee >(seqtk seq -1 - | pigz -c -p {threads} > {output.reads[0]}) | \
                        seqtk seq -2 - | pigz -c -p {threads} > {output.reads[1]}
                        ''')
            else:
                if reads_num == 1:
                    os.symlink(os.path.realpath(input[0]), output.reads[0])
                else:
                    shell('''cat {input} > {output.reads[0]}''')

        elif READS_FORMAT == "sra":
            if reads_num == 1:
                sra_file = os.path.basename(input[0])
                shell(
                    f'''
                     rm -rf {params.output_dir}/{sra_file}*
                     rm -rf {params.output_dir}.{sra_file}.temp

                     fasterq-dump \
                     --threads {threads} \
                     --split-3 \
                     --temp {params.output_dir}.{sra_file}.temp \
                     --outdir {params.output_dir} {input[0]} >{log} 2>&1

                     rm -rf {params.output_dir}.{sra_file}.temp
                     pigz --processes {threads} {params.output_dir}/{sra_file}_1.fastq
                     pigz --processes {threads} {params.output_dir}/{sra_file}_2.fastq
                     rm -rf {params.output_dir}/{sra_file}._*.fastq

                     mv {params.output_dir}/{sra_file}_1.fastq.gz {output.reads[0]}
                     mv {params.output_dir}/{sra_file}_2.fastq.gz {output.reads[1]}
                     ''')

            else:
                r1_list = []
                r2_list = []
                for sra in input:
                    sra_file = os.path.basename(sra)
                    r1_list.append(os.path.join(params.output_dir,
                                             sra_file + "_1.fastq.gz"))
                    r2_list.append(os.path.join(params.output_dir,
                                             sra_file + "_2.fastq.gz"))
                    shell(
                        f'''
                        rm -rf {params.output_dir}/{sra_file}*
                        rm -rf {params.output_dir}.{sra_file}.temp

                        fasterq-dump \
                        --threads {threads} \
                        --split-3 \
                        --temp {params.output_dir}.{sra_file}.temp \
                        --outdir {params.output_dir} {sra} >>{log} 2>&1

                        rm -rf {params.output_dir}.{sra_file}.temp
                        pigz --processes {threads} {params.output_dir}/{sra_file}_1.fastq
                        pigz --processes {threads} {params.output_dir}/{sra_file}_2.fastq
                        rm -rf {params.output_dir}/{sra_file}._*.fastq
                        ''')

                r1_str = " ".join(r1_list)
                r2_str = " ".join(r2_list)
                shell('''cat %s > %s''' % (r1_str, output.reads[0]))
                shell('''cat %s > %s''' % (r2_str, output.reads[1]))
                shell('''rm -rf %s''' % r1_str)
                shell('''rm -rf %s''' % r2_str)


rule prepare_short_reads_all:
    input:
        expand(os.path.join(
            config["output"]["raw"],
            "short_reads/{sample}/{sample}.raw{read}.fq.gz"),
            read=short_reads_suffix(),
            sample=SAMPLES.index.unique())


def get_reads(wildcards, step):
    read = short_reads_suffix()

    short_reads = expand(os.path.join(
        config["output"][step],
        "short_reads/{sample}/{sample}.{step}{read}.fq.gz"),
                         step=step,
                         read=read,
                         sample=wildcards.sample)
    return short_reads


if config["params"]["raw"]["fastqc"]["do"]:
    rule raw_fastqc:
        input:
            lambda wildcards: get_reads(wildcards, "raw", False)
        output:
            directory(os.path.join(
                config["output"]["raw"],
                "fastqc/{sample}.fastqc.out"))
        threads:
            config["params"]["raw"]["threads"]
        log:
            os.path.join(config["output"]["raw"], "logs/{sample}.fastqc.log")
        shell:
            '''
            mkdir -p {output}
            fastqc \
            --outdir {output} \
            --threads {threads} \
            --format fastq \
            {input} \
            2> {log}
            '''


    rule raw_fastqc_multiqc:
        input:
            expand(os.path.join(
                config["output"]["raw"],
                "fastqc/{sample}.fastqc.out"),
                   sample=SAMPLES.index.unique())
        output:
            html = os.path.join(
                config["output"]["raw"],
                "report/fastqc_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["raw"],
                "report/fastqc_multiqc_report_data"))
        params:
            output_dir = os.path.join(config["output"]["raw"],
                                      "report")
        log:
            os.path.join(config["output"]["raw"], "logs/multiqc_fastqc.log")
        shell:
            '''
            multiqc \
            --outdir {params.output_dir} \
            --title fastqc \
            --module fastqc {input} \
            2> {log}
            '''


    rule raw_fastqc_all:
        input:
            expand([
                os.path.join(
                    config["output"]["raw"],
                    "fastqc/{sample}.fastqc.out"),
                os.path.join(
                    config["output"]["raw"],
                    "report/fastqc_multiqc_report.html"),
                os.path.join(
                    config["output"]["raw"],
                    "report/fastqc_multiqc_report_data")],
                   sample=SAMPLES.index.unique())

else:
    rule raw_fastqc_all:
        input:


if config["params"]["qcreport"]["do"]:
    rule raw_report:
        input:
            lambda wildcards: get_reads(wildcards, "raw")
        output:
            os.path.join(config["output"]["raw"],
                         "report/stats/{sample}_raw_stats.tsv")
        params:
            sample_id = "{sample}",
            fq_encoding = config["params"]["fq_encoding"]
        log:
            os.path.join(config["output"]["raw"],
                         "logs/{sample}.seqkit.log")
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
                rnapi.change(output[0], params.sample_id, "raw",
                              "pe", ["fq1", "fq2"])
            else:
                rnapi.change(output[0], params.sample_id, "raw",
                              "se", ["fq1"])


    rule raw_report_merge:
        input:
            expand(
                os.path.join(config["output"]["raw"],
                             "report/stats/{sample}_raw_stats.tsv"),
                sample=SAMPLES.index.unique())
        output:
            os.path.join(config["output"]["qcreport"], "raw_stats.tsv")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            rnapi.merge(input, rnapi.parse, threads, output=output[0])


    rule raw_report_all:
        input:
            os.path.join(config["output"]["qcreport"], "raw_stats.tsv")

else:
    rule raw_report_all:
        input:


rule raw_all:
    input:
        #rules.prepare_reads_all.input,
        rules.raw_fastqc_all.input,
        rules.raw_report_all.input
