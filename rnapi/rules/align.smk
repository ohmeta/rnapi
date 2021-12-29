def get_clean_reads(wildcards):
    if TRIMMING_DO:
        return get_reads(wildcards, "trimming")
    else:
        return get_reads(wildcards, "raw")


rule align_star:
    input:
        reads = get_clean_reads,
        gtf = config["reference"]["gtf"],
        index = config["reference"]["index_star"]
    output:
        bam = os.path.join(config["output"]["align"], "star/{sample}/Aligned.sortedByCoord.out.bam"),
        tab = os.path.join(config["output"]["align"], "star/{sample}/ReadsPerGene.out.tab")
    params:
        outprefix = os.path.join(config["output"]["align"], "star/{sample}"),
        quant_mode = config["params"]["align"]["star"]["quant_mode"]
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/align_star/{sample}.log")
    run:
        shell(
            '''
            STAR \
            --quantMode {params.quant_mode} \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --genomeDir {input.index} \
            --readFilesIn {input.reads} \
            --outFileNamePrefix {params.outprefix}/ \
            --outSAMtype BAM SortedByCoordinate \
            --outStd Log \
            {log}
            ''')


if config["params"]["align"]["star"]["do"]:
    rule align_genome_star_all:
        input:
            expand([
                os.path.join(config["output"]["align"],
                             "star/{sample}/Aligned.sortedByCoord.out.bam"),
                os.path.join(config["output"]["align"],
                             "star/{sample}/ReadsPerGene.out.tab")],
                   sample=SAMPLES.index.unique())
else:
    rule align_genome_star_all:
        input:


rule align_all:
    input:
        rules.align_genome_star_all.input
