rule quantify_genome_star:
    input:
        expand(os.path.join(config["output"]["align"],
                            "star/genome/{sample}/ReadsPerGene.out.tab"),
               sample=SAMPLES.index.unique())
    output:
        counts = os.path.join(config["output"]["quantify"], "star_genome_counts.tsv")
    threads:
        config["params"]["quantify"]["threads"]
    log:
        os.path.join(config["output"]["quantify"], "logs/quantify_genome_star.log")
    params:
        strandedness = config["params"]["strandedness"]
    run:
        import numpy as np

        rnapi.parse_gene_tab_init(params.strandedness)
        count_df = rnapi.merge_cols(input, rnapi.parse_gene_tab, threads).fillna(0)
        count_df = count_df[count_df.apply(np.sum, axis=1)>0]
        count_df.reset_index().to_csv(output[0], sep="\t", index=False)


if config["params"]["align"]["star"]["do"]:
    if config["params"]["align"]["star"]["quant_mode"]["GeneCounts"]:
        rule quantify_genome_star_all:
            input:
                os.path.join(config["output"]["quantify"], "star_genome_counts.tsv")
    else:
        rule quantify_genome_star_all:
            input:
else:
    rule quantify_genome_star_all:
        input:


rule quantify_all:
    input:
        rules.quantify_genome_star_all.input
