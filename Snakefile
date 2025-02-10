configfile: "config.yaml"

# Include modules
include: "workflow/rules/preprocess.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/quant.smk"
#include: "workflow/rules/analysis.smk"

rule all:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
#        expand("05.Normalized/counts_normalized.tsv", **config),
#        expand("06.DEA/{comparison}_DE_results.tsv", comparison=config['diff_comparisons']),
#        expand("07.Visualization/venn_diagram.pdf", **config),
#        expand("08.Enrichment/GO_results.tsv", **config),
#        expand("08.Enrichment/KEGG_results.tsv", **config)

