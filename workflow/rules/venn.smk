rule generate_venn:
    input:
        de_files = expand("06.DGE_DESeq2/{comp}_DE_results.tsv", 
                        comp=config['diff_comparisons'].keys())
    output:
        "07.Visualization/venn_diagram.pdf"
    params:
        gene_sets = config['venn_genes']
    script:
        "../scripts/venn.R"