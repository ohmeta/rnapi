def get_clean_reads(wildcards):
    if TRIMMING_DO:
        return get_reads(wildcards, "trimming", False)
    else:
        return get_reads(wildcards, "raw", False)


rule index_star:
    input:
        dna = config["reference"]["dna"],
        gtf = config["reference"]["gtf"]
    output:
        directory(config["reference"]["index_star"])
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/star_index_genome.log")
    shell:
        '''
        mkdir -p {output}
        pigz -dk {input.dna} > {output}/genome.fasta
        pigz -dk {input.gtf} > {output}/genome.gtf

        STAR \
        --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output} \
        --genomeFastaFiles {output/genome.dna} \
        --sjdbGTFfile {output/genome.gtf} \
        --sjdbOverhang 100 \
        {log}

        rm -rf {output}/genome.fasta
        rm -rf {output}/genome.gtf
        '''


rule align_star:
    input:
        reads = get_clean_reads,
        gtf = config["reference"]["gtf"],
        index = config["reference"]["index_star"]
    output:
        bam = os.path.join(config["output"]["align"], "star/{sample}/Aligned.sortedByCoord.out.bam"),
        tab = os.path.join(config["output"]["align"], "star/{sample}/ReadsPerGene.out.tab")
    params:
        outprefix = os.path.join(config["output"]["align"], "star/{sample}")
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/align_star/{sample}.log")
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
            {log}
            ''')


if config["params"]["align"]["star"]["do"]:
    rule align_star_all:
        input:
            expand([
                os.path.join(config["output"]["align"],
                             "star/{sample}/Aligned.sortedByCoord.out.bam"),
                os.path.join(config["output"]["align"],
                             "star/{sample}/ReadsPerGene.out.tab")],
                   sample=SAMPLES.index.unique())
else:
    rule align_star_all:
        input:


rule align_all:
    rules.align_star_all.input
