rule quantify_gene_star:
    input:
        expand(os.path.join(config["output"]["align"],
                            "star/genome/{sample}/ReadsPerGene.out.tab"),
               sample=SAMPLES.index.unique())
    output:
        counts = os.path.join(config["output"]["quantify"], "star_gene_counts.tsv")
    threads:
        config["params"]["quantify"]["threads"]
    log:
        os.path.join(config["output"]["quantify"], "logs/quantify_gene_star.log")
    params:
        strandedness = config["params"]["strandedness"]
    run:
        import numpy as np

        rnapi.parse_gene_tab_init(params.strandedness)
        count_df = rnapi.merge_cols(input, rnapi.parse_gene_tab, threads).fillna(0)
        count_df = count_df[count_df.apply(np.sum, axis=1)>0]
        count_df.reset_index().to_csv(output[0], sep="\t", index=False)


rule quantify_transcript_star:
    input:
        bam = os.path.join(
            config["output"]["align"],
            "star/transcriptome/{sample}/Aligned.toTranscriptome.out.bam"),
        index = config["reference"]["index_rsem"]
    output:
        directory(os.path.join(config["output"]["quantify"], "star_transcript_counts/{sample}"))
    threads:
        config["params"]["quantify"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/quantify_transcript_star/{sample}.log")
    shell:
        '''
        rsem-calculate-expression \
        --bam --no-bam-output \
        -p {threads} --paired-end --forward-prob 0 \
        {input.bam} {input.index} {output} \
        > {log} 2>&1
        '''


if config["params"]["align"]["star"]["do"]:
    if config["params"]["align"]["star"]["quant_mode"]["GeneCounts"]:
        rule quantify_gene_star_all:
            input:
                os.path.join(config["output"]["quantify"], "star_gene_counts.tsv")
    else:
        rule quantify_gene_star_all:
            input:

    if config["params"]["align"]["star"]["quant_mode"]["TranscriptomeSAM"]:
        rule quantify_transcript_star_all:
            input:
                expand(os.path.join(config["output"]["quantify"], "star_transcript_counts/{sample}"),
                       sample=SAMPLES.index.unique())
    else:
        rule quantify_transcript_star_all:
            input:
else:
    rule quantify_gene_star_all:
        input:

    rule quantify_transcript_star_all:
        input:


rule quantify_all:
    input:
        rules.quantify_gene_star_all.input,
        rules.quantify_transcript_star_all.input
