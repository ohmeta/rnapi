rule index_star:
    input:
        dna = config["reference"]["dna"],
        gtf = config["reference"]["gtf"]
    output:
        expand(os.path.join(config["reference"]["index_star"], "{file}"),
               file=["Genome",  "SA", "SAindex"])
    params:
        index = config["reference"]["index_star"]
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/index_star.log")
    shell:
        '''
        mkdir -p {params.index}
        pigz -dkc {input.dna} > {params.index}/genome.fasta

        STAR \
        --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --genomeFastaFiles {params.index}/genome.fasta \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 100 \
        > {log} 2>&1

        rm -rf {params.index}/genome.fasta
        '''


rule index_rsem:
    input:
        dna = config["reference"]["dna"],
        gtf = config["reference"]["gtf"]
    output:
        expand(config["reference"]["index_rsem"] + ".{suffix}",
               suffix=["chrlist", "grp", "idx.fa", "n2g.idx.fa",
                       "seq", "ti", "transcripts.fa"])
    params:
        outprefix = config["reference"]["index_rsem"]
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/index_rsem.log")
    run:
        outdir = os.path.dirname(params.outprefix)

        shell(
            '''
            mkdir -p {outdir}
            pigz -dkc {input.dna} > {params.outprefix}.fasta

            rsem-prepare-reference \
            --gtf {input.gtf} \
            --num-threads {threads} \
            {params.outprefix}.fasta {params.outprefix} \
            > {log} 2>&1

            rm -rf {params.outprefix}.fasta
            ''')


rule index_salmon:
    input:
        dna = config["reference"]["dna_gencode"],
        cdna = config["reference"]["cdna_gencode"]
    output:
        expand(os.path.join(config["reference"]["index_salmon"], "{file}"),
               file=["complete_ref_lens.bin", "ctable.bin", "ctg_offsets.bin",
                     "mphf.bin", "pos.bin", "rank.bin", "refAccumLengths.bin",
                     "reflengths.bin", "refseq.bin", "seq.bin"])
    params:
        index = config["reference"]["index_salmon"],
        kmer_len = config["params"]["quantify"]["salmon"]["kmer_len"],
        index_add_genome = config["params"]["quantify"]["salmon"]["index_add_genome"]
    threads:
        config["params"]["align"]["threads"]
    log:
        os.path.join(config["output"]["align"], "logs/index_salmon.log")
    run:
        shell('''mkdir -p {params.index}''')

        if not params.index_add_genome:
            shell(
                '''
                salmon index \
                --transcripts {input.cdna} \
                --index {params.index} \
                --kmerLen {params.kmer_len} \
                --threads {threads} \
                > {log} 2>&1
                ''')
        else:
            shell(
                '''
                zcat {input.dna} | grep "^>" | awk -F'[> ]' '{{print $2}}' > {params.index}/decoys.txt
                zcat {input.cdna} {input.dna} > {params.index}/gentrome.fa.gz

                salmon index \
                --transcripts {params.index}/gentrome.fa.gz \
                --decoys {params.index}/decoys.txt \
                --index {params.index} \
                --kmerLen {params.kmer_len} \
                --threads {threads} \
                > {log} 2>&1
                ''')
