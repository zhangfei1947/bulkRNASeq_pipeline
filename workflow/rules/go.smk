rule go_enrichment:
    input:
        diff = "06.Diff_Expression/{comparison}.FC2.padj0.05.deseq2_results.tsv"
    output:
        outfile = "08.Enrichment/{comparison}.go_enrichment.tsv",
        pdf = "08.Enrichment/{comparison}.go_enrichment.pdf",
        png = "08.Enrichment/{comparison}.go_enrichment.png"
    log:
        "logs/go/{comparison}.log"
    script:
        "../scripts/go_enrich.R"
