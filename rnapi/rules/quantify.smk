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


# RSEM has the ability to produce both gene and isoform-level expression estimates.
# However, accurate isoform level expression is typically much more challenging than gene-level estimation,
# and isoform-level estimates are far noisier.
#
# These are tab-delimited files and contain expression estimates for each
# isoform ("transcript_id") or gene ("gene_id") as "expected_count", and also
# as TPM (Transcripts Per Million) and FPKM (Fragments Per Kilobase of transcript per Million mapped reads) numbers.

rule quantify_transcript_star:
    input:
        bam = os.path.join(
            config["output"]["align"],
            "star/transcriptome/{sample}/Aligned.toTranscriptome.out.bam"),
        index = expand(config["reference"]["index_rsem"] + ".{suffix}",
                       suffix=["chrlist", "grp", "idx.fa", "n2g.idx.fa",
                               "seq", "ti", "transcripts.fa"])
    output:
        genes = os.path.join(config["output"]["quantify"],
                             "star_transcript_counts/{sample}/{sample}.genes.results"),
        isoforms = os.path.join(config["output"]["quantify"],
                                "star_transcript_counts/{sample}/{sample}.isoforms.results"),
        stats_cnt = os.path.join(config["output"]["quantify"],
                                 "star_transcript_counts/{sample}/{sample}.stat/{sample}.cnt"),
        stats_model = os.path.join(config["output"]["quantify"],
                                   "star_transcript_counts/{sample}/{sample}.stat/{sample}.model"),
        stats_theta = os.path.join(config["output"]["quantify"],
                                   "star_transcript_counts/{sample}/{sample}.stat/{sample}.theta")
    params:
        strandedness = config["params"]["strandedness"],
        index = config["reference"]["index_rsem"],
        outprefix = os.path.join(config["output"]["quantify"],
                                 "star_transcript_counts/{sample}/{sample}")
    threads:
        config["params"]["quantify"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/quantify_transcript_star/{sample}.log")
    run:
        import pandas as pd

        forward_prob  = 0
        if pd.isnull(params.strandedness) or \
           params.strandedness == "none" or \
           params.strandedness == "":
            forward_prob = 0.5  # non stranded protocol
        elif strandedness == "forward":
            forward_prob = 1 # for ligation-stranded protocols
        elif strandedness == "reverse":
            forward_prob = 0 # for dUTP libraries
        else:
            sys.exit("strandedness is not right")

        shell(
            '''
            rsem-calculate-expression \
            --bam --no-bam-output \
            -p {threads} \
            --paired-end \
            --forward-prob {forward_prob} \
            {input.bam} {params.index} {params.outprefix} \
            > {log} 2>&1
            ''')


rule quantify_transcript_star_merge:
    input:
        genes = expand(os.path.join(config["output"]["quantify"],
                             "star_transcript_counts/{sample}/{sample}.genes.results"),
                       sample=SAMPLES.index.unique()),
        transcripts = expand(os.path.join(config["output"]["quantify"],
                                          "star_transcript_counts/{sample}/{sample}.isoforms.results"),
                             sample=SAMPLES.index.unique())
    output:
        gene_tpm = os.path.join(config["output"]["quantify"], "star_gene_counts_TPM.tsv"),
        gene_fpkm = os.path.join(config["output"]["quantify"], "star_gene_counts_FPKM.tsv"),
        transcript_tpm = os.path.join(config["output"]["quantify"], "star_transcript_counts_TPM.tsv"),
        transcript_fpkm = os.path.join(config["output"]["quantify"], "star_transcript_counts_FPKM.tsv")
    threads:
        config["params"]["quantify"]["threads"]
    run:
        import numpy as np

        def save_df(input_list, threads, func, outf):
            df = rnapi.merge_cols(input_list, func, threads).fillna(0)
            df = df[df.apply(np.sum, axis=1)>0]
            df.reset_index().to_csv(outf, sep="\t", index=False)

        save_df(input.genes, threads, rnapi.parse_rsem_gene_TPM, output.gene_tpm)
        save_df(input.genes, threads, rnapi.parse_rsem_gene_FPKM, output.gene_fpkm)
        save_df(input.transcripts, threads, rnapi.parse_rsem_transcript_TPM, output.transcript_tpm)
        save_df(input.transcripts, threads, rnapi.parse_rsem_transcript_FPKM, output.transcript_fpkm)


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
                os.path.join(config["output"]["quantify"], "star_gene_counts_TPM.tsv"),
                os.path.join(config["output"]["quantify"], "star_gene_counts_FPKM.tsv"),
                os.path.join(config["output"]["quantify"], "star_transcript_counts_TPM.tsv"),
                os.path.join(config["output"]["quantify"], "star_transcript_counts_FPKM.tsv")
    else:
        rule quantify_transcript_star_all:
            input:
else:
    rule quantify_gene_star_all:
        input:

    rule quantify_transcript_star_all:
        input:



rule quantify_salmon:
    input:
        reads = get_clean_reads,
        gtf = config["reference"]["gtf"],
        index = expand(os.path.join(config["reference"]["index_salmon"], "{file}"),
                       file=["complete_ref_lens.bin", "ctable.bin", "ctg_offsets.bin",
                             "mphf.bin", "pos.bin", "rank.bin", "refAccumLengths.bin",
                             "reflengths.bin", "refseq.bin", "seq.bin"])
    output:
        gene_qf = os.path.join(config["output"]["quantify"], "salmon/{sample}/quant.genes.sf"),
        transcript_qf = os.path.join(config["output"]["quantify"], "salmon/{sample}/quant.sf")
    params:
        index = config["reference"]["index_salmon"],
        outdir = os.path.join(config["output"]["quantify"], "salmon/{sample}"),
        lib_type = config["params"]["quantify"]["salmon"]["lib_type"],
        extra = config["params"]["quantify"]["salmon"]["extra"]
    threads:
        config["params"]["quantify"]["threads"]
    shell:
        '''
        salmon quant \
        --index {params.index} \
        --geneMap {input.gtf} \
        --libType {params.lib_type} \
        -1 {input.reads[0]} \
        -2 {input.reads[1]} \
        --output {params.outdir} \
        --threads {threads} \
        {params.extra}
        '''


rule quantify_salmon_merge:
    input:
        genes = expand(os.path.join(config["output"]["quantify"], "salmon/{sample}/quant.genes.sf"),
                       sample=SAMPLES.index.unique()),
        transcripts = expand(os.path.join(config["output"]["quantify"], "salmon/{sample}/quant.sf"),
                             sample=SAMPLES.index.unique())
    output:
        gene_count = os.path.join(config["output"]["quantify"], "salmon_gene_counts.tsv"),
        gene_tpm = os.path.join(config["output"]["quantify"], "salmon_gene_counts_TPM.tsv"),
        transcript_count = os.path.join(config["output"]["quantify"], "salmon_transcript_counts.tsv"),
        transcript_tpm = os.path.join(config["output"]["quantify"], "salmon_transcript_counts_tpm.tsv")
    threads:
        config["params"]["quantify"]["threads"]
    run:
        import numpy as np

        def quant_merger(input_list, func):
            df = rnapi.merge_cols(input_list, func, threads).fillna(0)
            df = df[df.apply(np.sum, axis=1)>0]
            return df.reset_index()

        quant_merger(input.genes, rnapi.parse_salmon_TPM)\
            .rename(columns={"Name": "gene_id"})\
            .to_csv(output.gene_tpm, sep="\t", index=False)

        quant_merger(input.genes, rnapi.parse_salmon_count)\
            .rename(columns={"Name": "gene_id"})\
            .to_csv(output.gene_count, sep="\t", index=False)

        quant_merger(input.transcripts, rnapi.parse_salmon_TPM)\
            .rename(columns={"Name": "transcript_id"})\
            .to_csv(output.transcripts_tpm, sep="\t", index=False)

        quant_merger(input.transcripts, rnapi.parse_salmon_count)\
            .rename(columns={"Name": "transcript_id"})\
            .to_csv(output.transcripts_count, sep="\t", index=False)


if config["params"]["quantify"]["salmon"]["do"]:
    rule quantify_salmon_all:
        input:
            os.path.join(config["output"]["quantify"], "salmon_gene_counts.tsv"),
            os.path.join(config["output"]["quantify"], "salmon_gene_counts_TPM.tsv"),
            os.path.join(config["output"]["quantify"], "salmon_transcript_counts.tsv"),
            os.path.join(config["output"]["quantify"], "salmon_transcript_counts_tpm.tsv")

else:
    rule quantify_salmon_all:
        input:


rule quantify_all:
    input:
        rules.quantify_gene_star_all.input,
        rules.quantify_transcript_star_all.input,
        rules.quantify_salmon_all.input
