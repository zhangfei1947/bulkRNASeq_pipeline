rule differential_expression:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
    output:
        "06.DGE_DESeq2/{comparison}_results.tsv"
    params:
        contrast = lambda wildcards: config['diff_comparisons'][wildcards.comparison]
    log:
        "logs/diff/{comparison}.log"
    script:
        "../scripts/diff.R"