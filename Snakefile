envmodules: True

#configfile: "config.yaml"

# Include modules
# softlink rawdata 
include: "workflow/rules/preprocess.smk"

# rawdata quality control
include: "workflow/rules/qc.smk"

# mapping to reference genome
include: "workflow/rules/align.smk"

# gene expression quantification
include: "workflow/rules/quant.smk"

# row readcount normalization, sample correlation, PCA
#include: "workflow/rules/normalize.smk"

# differential gene expression
#include: "workflow/rules/dge.smk"

# plot venn diagram
#include: "workflow/rules/venn.smk"

# functional analysis
#include: "workflow/rules/func.smk"


rule all:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
#        expand("05.Normalized/counts_normalized.tsv", **config),
#        expand("06.DEA/{comparison}_DE_results.tsv", comparison=config['diff_comparisons']),
#        expand("07.Visualization/venn_diagram.pdf", **config),
#        expand("08.Enrichment/GO_results.tsv", **config),
#        expand("08.Enrichment/KEGG_results.tsv", **config)

