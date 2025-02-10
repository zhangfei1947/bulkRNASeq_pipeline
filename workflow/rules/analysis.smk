rule normalize_counts:
    input:
        "04.Quantification/counts_raw.tsv"
    output:
        normalized = "05.Normalized/counts_normalized.tsv",
        fpkm = "05.Normalized/fpkm_matrix.tsv"
    script:
        "../scripts/norm.R"

rule differential_expression:
    input:
        "04.Quantification/counts_raw.tsv"
    output:
        "06.DEA/{comparison}_DE_results.tsv"
    params:
        contrast = lambda wildcards: config['diff_comparisons'][wildcards.comparison]
    log:
        "logs/dea/{comparison}.log"
    script:
        "../scripts/dea.R"

rule generate_venn:
    input:
        de_files = expand("06.DEA/{comp}_DE_results.tsv", 
                        comp=config['diff_comparisons'].keys())
    output:
        "07.Visualization/venn_diagram.pdf"
    params:
        gene_sets = config['venn_genes']
    script:
        "../scripts/venn.R"

rule go_enrichment:
    input:
        de_files = expand("06.DEA/{comp}_DE_results.tsv",
                        comp=config['diff_comparisons'].keys())
    output:
        "08.Enrichment/GO_results.tsv",
        "08.Enrichment/GO_plots.pdf"
    params:
        databases = config['enrichment']['GO']['databases'],
        fdr = config['enrichment']['GO']['fdr_cutoff']
    script:
        "../scripts/go_enrich.R"

rule kegg_enrichment:
    input:
        de_files = expand("06.DEA/{comp}_DE_results.tsv",
                        comp=config['diff_comparisons'].keys())
    output:
        "08.Enrichment/KEGG_results.tsv",
        "08.Enrichment/KEGG_plots.pdf"
    params:
        org_db = config['enrichment']['KEGG']['org_db'],
        fdr = config['enrichment']['KEGG']['fdr_cutoff']
    script:
        "../scripts/kegg_enrich.R"

rule analysis_summary:
    input:
        normalized = rules.normalize_counts.output,
        dea = expand("06.DEA/{comp}_DE_results.tsv", 
                   comp=config['diff_comparisons'].keys()),
        venn = rules.generate_venn.output,
        go = rules.go_enrichment.output,
        kegg = rules.kegg_enrichment.output
    output:
        "09.Reports/final_summary.md"
    run:
        template = """
        # RNA-Seq Analysis Summary

        ## Normalization
        - Normalized counts: {normalized}
        - FPKM matrix: {fpkm}

        ## Differential Expression
        {dea_list}

        ## Enrichment Results
        - Significant GO terms: {go_count}
        - Significant KEGG pathways: {kegg_count}
        """
        dea_list = "\n".join([f"- {comp}" for comp in config['diff_comparisons']])
        
        with open(output[0], 'w') as f:
            f.write(template.format(
                normalized=input.normalized[0],
                fpkm=input.normalized[1],
                dea_list=dea_list,
                go_count=sum(1 for line in open(input.go[0]) if line.startswith("GO")),
                kegg_count=sum(1 for line in open(input.kegg[0]) if not line.startswith("#"))
            ))

