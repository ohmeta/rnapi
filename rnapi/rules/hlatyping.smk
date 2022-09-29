rule hlatyping_arcashla_reference:
    output:
        os.path.join(config["output"]["hlatyping"], "config/arcasHLA_reference_done")
    log:
        os.path.join(config["output"]["hlatyping"], "logs/arcasHLA_reference.log")
    params:
        version = config["params"]["hlatyping"]["arcashla"]["IMGTHLA_version"]
    conda:
        config["envs"]["arcashla"]
    shell:
        '''
        if [ ! -e {output} ];
        then
            if [ "{params.version}" == "latest" ];
            then
                git lfs install >{log} 2>&1
                arcasHLA reference --update -v >>{log} 2>&1
                touch {output}
            else
                git lfs install >{log} 2>&1
                arcasHLA reference --version {params.version} -v >>{log} 2>&1
                touch {output}
            fi
        fi
        '''


rule hlatyping_arcashla_extract:
    input:
        done = os.path.join(config["output"]["hlatyping"], "config/arcasHLA_reference_done"),
        bam = os.path.join(config["output"]["align"],
                           "star/reads/{sample}/Aligned.sortedByCoord.out.bam")
    output:
        r1 = os.path.join(config["output"]["hlatyping"], "reads/{sample}/{sample}.extracted.1.fq.gz"),
        r2 = os.path.join(config["output"]["hlatyping"], "reads/{sample}/{sample}.extracted.2.fq.gz")
    log:
        os.path.join(config["output"]["hlatyping"],
                     "logs/arcashla_extract/{sample}.archashla_extract.log")
    benchmark:
        os.path.join(config["output"]["hlatyping"],
                     "benchmark/arcashla_extract/{sample}.archashla_extract.benchmark.txt")
    conda:
        config["envs"]["arcashla"]
    params:
        outdir = os.path.join(config["output"]["hlatyping"], "reads/{sample}"),
        unmapped = "--unmapped" if config["params"]["hlatyping"]["arcashla"]["unmapped"] else ""
    threads:
        config["params"]["hlatyping"]["threads"]
    shell:
        '''
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        arcasHLA extract \
        --threads {threads} \
        {params.unmapped} \
        --outdir {params.outdir} \
        --verbose \
        --log {log} \
        {input.bam}
        '''


rule hlatyping_arcashla_genotype:
    input:
        done = os.path.join(config["output"]["hlatyping"], "config/arcasHLA_reference_done"),
        r1 = os.path.join(config["output"]["hlatyping"], "reads/{sample}/{sample}.extracted.1.fq.gz"),
        r2 = os.path.join(config["output"]["hlatyping"], "reads/{sample}/{sample}.extracted.2.fq.gz")
    output:
        json = os.path.join(config["output"]["hlatyping"], "genotype/{sample}/{sample}.genotype.json")
    log:
        os.path.join(config["output"]["hlatyping"],
                     "logs/arcashla_genotype/{sample}.archashla_genotype.log")
    benchmark:
        os.path.join(config["output"]["hlatyping"],
                     "benchmark/arcashla_genotype/{sample}.archashla_genotype.benchmark.txt")
    conda:
        config["envs"]["arcashla"]
    params:
        genes = ",".join(config["params"]["hlatyping"]["arcashla"]["genes"]),
        outdir = os.path.join(config["output"]["hlatyping"], "genotype/{sample}")
    threads:
        config["params"]["hlatyping"]["threads"]
    shell:
        '''
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        arcasHLA genotype \
        {input} \
        --genes {params.genes} \
        --outdir {params.outdir} \
        --threads {threads} \
        --verbose \
        --log {log}
        '''


if config["params"]["hlatyping"]["arcashla"]["do"]:
    rule hlatyping_arcashla_all:
        input:
            expand(
                os.path.join(config["output"]["hlatyping"], "genotype/{sample}/{sample}.genotype.json"),
                sample=SAMPLES.index.unique())
else:   
    rule hlatyping_arcashla_all:
        input:


rule hlatyping_all:
    input:
        rules.hlatyping_arcashla_all.input
            

localrules:
    hlatyping_arcashla_reference