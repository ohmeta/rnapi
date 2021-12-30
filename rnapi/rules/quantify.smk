rule quantify_genome_star:
    input:
        expand(
            os.path.join(config["output"]["align"],
                         "star/genome/{sample}/ReadsPerGene.out.tab"),
            sample=SAMPLES.index.unique())
    output:
        counts = os.path.join(config["output"]["quantify"], "star_genome_counts.tsv")
    threads:
        config["output"]["quantify"]["threads"]
    log:
        os.path.join(config["output"]["quantify"], "logs/quantify_genome_star.log")
    params:
        strandedness = config["params"]["strandedness"]
    run:
        global STAR_GENE_TAB_COLUMN__ = 3

        if pd.isnull(params.strandedness) or \
           params.strandedness == "none" or \
           params.strandedness == "":
            STAR_GENE_TAB_COLUMN__  = 1  # non stranded protocol
        elif strandedness == "yes":
            STAR_GENE_TAB_COLUMN__ = 2  # 3rd column
        elif strandedness == "reverse":
            STAR_GENE_TAB_COLUMN__ = 3  # 4th column, usually for Illumina truseq
        else:
            sys.exit("strandedness is not right")

        def parse_gene_tab(genef):
            sample_id = os.path.basename(os.path.dirname(genef))
            if os.path.exists(genef):
                try:
                    df = pd.read_csv(genef, index_col=0, usecols=[0, column], names=[sample_id], skiprows=4, sep="\t")
                except pd.errors.EmptyDataError:
                    print(f"{genef} is empty, please check")
                    return None

                if not df.empty:
                    return df
                else:
                    print(f"{genef} is empty, please check")
                    return None
          else:
              sys.exit(f"{genef} is not exists")

        count_df = rnapi.merge_cols(input, parse_gene_tab, threads)

        count_df.index.name = "gene"
        count_df = count_df.fillna(0)
        count_df = count_df[count_df.apply(np.sum, axis=1)>0]
        count_df.reset_index().to_csv(output[0], sep="\t")


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
