rule index_star:
    input:
        dna = config["reference"]["dna"],
        gtf = config["reference"]["gtf"]
    output:
        directory(config["reference"]["index_star"])
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/index_star.log")
    shell:
        '''
        mkdir -p {output}
        pigz -dkc {input.dna} > {output}/genome.fasta

        STAR \
        --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output} \
        --genomeFastaFiles {output}/genome.fasta \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 100 \
        {log}

        rm -rf {output}/genome.fasta
        '''


rule index_rsem:
    input:
        dna = config["reference"]["dna"],
        gtf = config["reference"]["gtf"]
    output:
        expand(config["reference"]["index_rsem"] + ".{suffix}",
               suffix=["chrlist", "grp", "idx.fa", "n2g.idx.fa", "seq", "ti", "transcripts.fa"])
    params:
        outprefix = config["reference"]["index_rsem"]
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/index_rsem.log")
    shell:
        '''
        mkdir -p {output}
        pigz -dkc {input.dna} > {params.outprefix}.fasta

        rsem-prepare-reference \
        --gtf {input.gtf} \
        --num-threads {threads} \
        {params.outprefix}.fasta {params.outprefix} \
        > {log} 2>&1

        rm -rf {params.outprefix}.fasta
        '''
