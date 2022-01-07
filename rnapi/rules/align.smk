# Reference
# Dobin, A. and Gingeras, T.R. 2015. Mapping RNA-seq reads with STAR. Curr. Protoc. Bioinform. 51:11.14.1-11.14.19. doi: 10.1002/0471250953.bi1114s51
#
# Protocol list
## 1. mapping RNA-seq reads to the reference genome
## 2. mapping RNA-seq reads with 2-pass procedure
## 3. mapping reads and generating unsorted and coordinate-sorted BAM files
## 4. generating signal files for visualization on genome browsers for stranded RNA-seq data
## 5. generating signal files for visualization on genome browsers for un-stranded RNA-seq data
## 6. mapping RNA-seq reads and generating chimeric alignments to detect fusion transcripts and circular RNA
## 7. mapping RNA-seq reads, generating output in transcriptomic coordinates and using RSEM to quantify expression of transcripts and genes
## 8. mapping RNA-seq reads and running cufflinks to assemble and quantify transcripts for stranded RNA-seq data
## 9. mapping RNA-seq reads and running cufflinks to assemble and quantify transcripts for un-stranded RNA-seq data


def get_clean_reads(wildcards):
    if DELRIBORNA_DO:
        return get_reads(wildcards, "delriborna")
    elif TRIMMING_DO:
        return get_reads(wildcards, "trimming")
    else:
        return get_reads(wildcards, "raw")


rule align_genome_star:
    input:
        reads = get_clean_reads,
        gtf = config["reference"]["gtf"],
        index = expand(os.path.join(config["reference"]["index_star"], "{file}"),
                       file=["Genome", "SA", "SAindex"])
    output:
        align_bam = os.path.join(config["output"]["align"], "star/genome/{sample}/Aligned.sortedByCoord.out.bam"),
        gene_tab = os.path.join(config["output"]["align"], "star/genome/{sample}/ReadsPerGene.out.tab")
    params:
        index = config["reference"]["index_star"],
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
            --genomeDir {params.index} \
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
        index = expand(os.path.join(config["reference"]["index_star"], "{file}"),
                       file=["Genome", "SA", "SAindex"])
    output:
        align_bam = os.path.join(config["output"]["align"],
                                 "star/transcriptome/{sample}/Aligned.sortedByCoord.out.bam"),
        trans_bam = os.path.join(config["output"]["align"],
                                 "star/transcriptome/{sample}/Aligned.toTranscriptome.out.bam")
    params:
        index = config["reference"]["index_star"],
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
            --genomeDir {params.index} \
            --readFilesIn {input.reads} \
            --outFileNamePrefix {params.outprefix}/ \
            --outSAMtype BAM SortedByCoordinate \
            --outStd Log \
            > {log} 2>&1
            ''')


if config["params"]["align"]["star"]["do"]:
    if config["params"]["align"]["star"]["quant_mode"]["GeneCounts"]:
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
