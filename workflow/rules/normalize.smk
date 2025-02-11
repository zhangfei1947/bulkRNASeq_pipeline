
rule normalize_counts:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
    output:
        normalized = "05.Normalization_DESeq2/counts_normalized.tsv",
        fpkm = "05.Normalization_DESeq2/fpkm_matrix.tsv"
    resources:
        mem_mb = 2000,
        runtime = 30,
        nodes = 1,
        ntasks_per_node = 1
    script:
        "../scripts/norm.R "