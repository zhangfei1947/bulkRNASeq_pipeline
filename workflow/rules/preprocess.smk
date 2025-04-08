from snakemake.utils import validate

#validate(config, schema="../schemas/config.schema.yaml")

def get_final_output():
    final_output = []
    # 01.raw data
    final_output.extend(expand("01.Rawdata/{sample}/{sample}_R1.fq.gz", sample=config["samples"]))
    final_output.extend(expand("01.Rawdata/{sample}/{sample}_R2.fq.gz", sample=config["samples"]))
    # 02.qc
    final_output.extend(expand("02.QC_fastp/{sample}/{sample}_R1.clean.fastq.gz", sample=config["samples"]))
    final_output.extend(expand("02.QC_fastp/{sample}/{sample}_R2.clean.fastq.gz", sample=config["samples"]))
    final_output.append("02.QC_fastp/duprate.boxplot.png")
    # 03.Alignment
    final_output.extend(expand("03.Alignment_hisat2/{sample}/{sample}.bam", sample=config["samples"]))
    final_output.extend(expand("03.Alignment_hisat2/{sample}/{sample}.summary", sample=config["samples"]))
    final_output.append("03.Alignment_hisat2/mapping.summary")
    final_output.append("03.Alignment_hisat2/mappingrate.boxplot.png")
    # 04.Quantification
    final_output.append("04.Quant_featureCounts/counts_raw.tsv")
    final_output.append("04.Quant_featureCounts/counts_raw.tsv.summary")
    final_output.append("04.Quant_featureCounts/counts_filter.tsv")
    final_output.append("04.Quant_featureCounts/fc.summary")
    final_output.append("04.Quant_featureCounts/assignment_stacked_barplot.png")
    # 05.Normalization
    final_output.append("05.Normalization_DESeq2/counts_normalized.tsv")
    final_output.append("05.Normalization_DESeq2/fpkm_matrix.tsv")
    if "corr" in config:
        final_output.extend(expand("05.Normalization_DESeq2/corr.heatmap.{corr_name}.png", corr_name=config["corr"]))
    if "pca_color" in config:
        final_output.extend(expand("05.Normalization_DESeq2/pca.plot.{color_scheme}.png", color_scheme=config["pca_color"]))
    # 06.DGE
    if "diff_comparisons" in config:
        final_output.extend(expand("06.Diff_Expression/{comparison}.deseq2_results.tsv", comparison=config['diff_comparisons']))
        final_output.extend(expand("06.Diff_Expression/{comparison}.FC2.padj0.05.deseq2_results.tsv", comparison=config['diff_comparisons']))
    # 07.Venn
    if "venn" in config:
        final_output.extend(expand("07.Venn/{group}.venn.png", group=config['venn']))
    # 08.Enrichment
    if "GO" in config:
        final_output.extend(expand("08.Enrichment/{comparison}.go_enrichment.tsv", comparison=config["GO"]))
        final_output.extend(expand("08.Enrichment/{comparison}.go_enrichment.pdf", comparison=config["GO"]))

    return final_output


localrules: make_links

rule make_links:
    input:
        r1 = lambda wildcards: config['samples'][wildcards.sample]['read1'],
        r2 = lambda wildcards: config['samples'][wildcards.sample]['read2']
    output:
        r1 = "01.Rawdata/{sample}/{sample}_R1.fq.gz",
        r2 = "01.Rawdata/{sample}/{sample}_R2.fq.gz"
    log:
        "logs/prep/{sample}.log"
    container: None
    shell:
        """
        ln -sf {input.r1} {output.r1} > {log} 2>&1
        ln -sf {input.r2} {output.r2} > {log} 2>&1
        """
