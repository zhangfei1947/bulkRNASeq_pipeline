rule go_enrichment:
    input:
        diff = "06.Diff_Expression/{comparison}.FC2.padj0.05.deseq2_results.tsv"
    output:
        "08.Enrichment/{comparison}.go_enrichment.tsv"
    envmodules:
        "GCC/12.2.0",
        "OpenMPI/4.1.4",
        "R/4.3.1"
    script:
        "../../scripts/go_enrich.R"
