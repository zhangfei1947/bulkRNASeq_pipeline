localrules: corr_heat, pca

rule normalize_counts:
    input:
        "04.Quant_featureCounts/counts_filter.tsv"
    output:
        normalized = "05.Normalization_DESeq2/counts_normalized.tsv",
        fpkm = "05.Normalization_DESeq2/fpkm_matrix.tsv"
    params:
        anno = config['genome']['geneloc'],
        sample = ",".join( config['samples'].keys() ),
        group = ",".join( sample['group'] for sample in config["samples"].values() ),
    script:
        "../scripts/norm.R"


rule corr_heat:
    input:
        norm_counts = "05.Normalization_DESeq2/counts_normalized.tsv"
    output:
        "05.Normalization_DESeq2/corr.heatmap.{corr_name}.png"
    params:
        target_groups = lambda wildcards: config["corr"][wildcards.corr_name].split(","),
        sample_mapping = lambda wildcards: {s: info["group"] for s, info in config["samples"].items()}
    script:
        "../scripts/corr.py"


rule pca:
    input:
        norm_counts = "05.Normalization_DESeq2/counts_normalized.tsv"
    output:
        "05.Normalization_DESeq2/pca.plot.{color_scheme}.png"
    params:
        color_mapping = lambda wildcards: {
            color_group: groups.split(",")
            for color_group, groups in config["pca_color"][wildcards.color_scheme].items()
        },

        sample_colors = lambda wildcards: {
            sample: next(
                (cg for cg, gs in config["pca_color"][wildcards.color_scheme].items() 
                 if info["group"] in gs.split(",")),
                "unknown" 
            )
            for sample, info in config["samples"].items()
        }
    script:
        "../scripts/pca.py"
