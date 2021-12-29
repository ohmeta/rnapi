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
        pigz -dk {input.dna} > {output}/genome.fasta
        pigz -dk {input.gtf} > {output}/genome.gtf

        STAR \
        --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output} \
        --genomeFastaFiles {output}/genome.dna \
        --sjdbGTFfile {output}/genome.gtf \
        --sjdbOverhang 100 \
        {log}

        rm -rf {output}/genome.fasta
        rm -rf {output}/genome.gtf
        '''
