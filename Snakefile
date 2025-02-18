#envmodules: True

#configfile: "config.yaml"

rule all:
    input:
        "02.QC_fastp/QC_summary_table.tsv",
        "02.QC_fastp/duprate.boxplot.png",
        "03.Alignment_hisat2/mapping.summary",
        "03.Alignment_hisat2/mappingrate.boxplot.png",
        "04.Quant_featureCounts/counts_raw.tsv",
        "04.Quant_featureCounts/counts_filter.tsv",
        "04.Quant_featureCounts/fc.summary",
        "04.Quant_featureCounts/assignment_stacked_barplot.png",
        "05.Normalization_DESeq2/counts_normalized.tsv",
        expand("05.Normalization_DESeq2/corr.heatmap.{corr_name}.png", corr_name=config["corr"].keys()),
        expand("05.Normalization_DESeq2/pca.plot.{color_scheme}.png", color_scheme=config["pca_color"].keys())
#        expand("05.Normalized/counts_normalized.tsv", **config),
#        expand("06.DEA/{comparison}_DE_results.tsv", comparison=config['diff_comparisons']),
#        expand("07.Visualization/venn_diagram.pdf", **config),
#        expand("08.Enrichment/GO_results.tsv", **config),
#        expand("08.Enrichment/KEGG_results.tsv", **config)



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
include: "workflow/rules/normalize.smk"

# differential gene expression
#include: "workflow/rules/dge.smk"

# plot venn diagram
#include: "workflow/rules/venn.smk"

# functional analysis
#include: "workflow/rules/func.smk"

