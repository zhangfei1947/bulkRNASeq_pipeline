rule normalize_counts:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
    output:
        normalized = "05.Normalization_DESeq2/counts_normalized.tsv",
        fpkm = "05.Normalization_DESeq2/fpkm_matrix.tsv"
    script:
        "../scripts/norm.R"