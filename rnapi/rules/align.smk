def get_clean_reads(wildcards):
    if TRIMMING_DO:
        return get_reads(wildcards, "trimming")
    else:
        return get_reads(wildcards, "raw")


rule align_genome_star:
    input:
        reads = get_clean_reads,
        gtf = config["reference"]["gtf"],
        index = config["reference"]["index_star"]
    output:
        align_bam = os.path.join(config["output"]["align"], "star/genome/{sample}/Aligned.sortedByCoord.out.bam"),
        gene_tab = os.path.join(config["output"]["align"], "star/genome/{sample}/ReadsPerGene.out.tab")
    params:
        outprefix = os.path.join(config["output"]["align"], "star/genome/{sample}")
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/align_genome_star/{sample}.log")
    run:
        shell(
            '''
            STAR \
            --quantMode GeneCounts \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --genomeDir {input.index} \
            --readFilesIn {input.reads} \
            --outFileNamePrefix {params.outprefix}/ \
            --outSAMtype BAM SortedByCoordinate \
            --outStd Log \
            > {log} 2>&1
            ''')


rule align_transcriptome_star:
    input:
        reads = get_clean_reads,
        gtf = config["reference"]["gtf"],
        index = config["reference"]["index_star"]
    output:
        align_bam = os.path.join(config["output"]["align"], "star/transcriptome/{sample}/Aligned.sortedByCoord.out.bam"),
        trans_bam = os.path.join(config["output"]["align"], "star/transcriptome/{sample}/Aligned.toTranscriptome.out.bam")
    params:
        outprefix = os.path.join(config["output"]["align"], "star/transcriptome/{sample}")
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/align_transcriptome_star/{sample}.log")
    run:
        shell(
            '''
            STAR \
            --quantMode TranscriptomeSAM \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --genomeDir {input.index} \
            --readFilesIn {input.reads} \
            --outFileNamePrefix {params.outprefix}/ \
            --outSAMtype BAM SortedByCoordinate \
            --outStd Log \
            > {log} 2>&1
            ''')


if config["params"]["align"]["star"]["do"]:
    if config["params"]["align"]["star"]["quant_mode"]["GeneCounts"]:
        print("GeneCounts: Yes")
        rule align_genome_star_all:
            input:
                expand([
                    os.path.join(config["output"]["align"],
                                 "star/genome/{sample}/Aligned.sortedByCoord.out.bam"),
                    os.path.join(config["output"]["align"],
                                 "star/genome/{sample}/ReadsPerGene.out.tab")],
                       sample=SAMPLES.index.unique())
    else:
        rule align_genome_star_all:
            input:


    if config["params"]["align"]["star"]["quant_mode"]["TranscriptomeSAM"]:
        print("TranscriptomeSAM: Yes")
        rule align_transcriptome_star_all:
            input:
                expand([
                    os.path.join(config["output"]["align"],
                                 "star/transcriptome/{sample}/Aligned.sortedByCoord.out.bam"),
                    os.path.join(config["output"]["align"],
                                 "star/transcriptome/{sample}/Aligned.toTranscriptome.out.bam")],
                       sample=SAMPLES.index.unique())
    else:
        rule align_transcriptome_star_all:
            input:

else:
    rule align_genome_star_all:
        input:


    rule align_transcriptome_star_all:
        input:


rule align_star_all:
    input:
        rules.align_genome_star_all.input,
        rules.align_transcriptome_star_all.input


rule align_all:
    input:
        rules.align_star_all.input
