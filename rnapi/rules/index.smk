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
        directory(config["reference"]["index_rsem"])
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/index_rsem.log")
    shell:
        '''
        mkdir -p {output}
        pigz -dkc {input.dna} > {output}/genome.fasta

        rsem-prepare-reference \
        --gtf {input.gtf} \
        {output}/genome.fasta {output} \
        > {log} 2>&1

        rm -rf {output}/genome.fasta
        '''
