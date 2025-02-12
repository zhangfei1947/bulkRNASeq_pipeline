
rule normalize_counts:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
    output:
        normalized = "05.Normalization_DESeq2/counts_normalized.tsv",
        fpkm = "05.Normalization_DESeq2/fpkm_matrix.tsv"
    resources:
        runtime = 30,
        nodes = 1,
        ntasks = 1,
        cpus_per_task= 1,
        mem_mb = 4000
    script:
        "../scripts/norm.R "