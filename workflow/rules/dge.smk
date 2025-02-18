

rule diff_analysis:
    input:
        counts = "04.Quant_featureCounts/counts_filter.tsv"
    output:
        results = "06.Diff_Expression/{comparison}/{comparison}.deseq2_results.tsv",
        sigresults = "06.Diff_Expression/{comparison}/{comparison}.deseq2_results.logFC1.padj0.05.tsv",
        vcpdf = "06.Diff_Expression/{comparison}/{comparison}.volcano.pdf",
        vcpng = "06.Diff_Expression/{comparison}/{comparison}.volcano.png"
    params:
        sample_info = lambda wildcards: {
            sample: info["group"]
            for sample, info in config["samples"].items()
            if info["group"] in [
                config["diff_comparisons"][wildcards.comparison]["numerator"],
                config["diff_comparisons"][wildcards.comparison]["denominator"]
            ]
        },
        comparison = lambda wildcards: config["diff_comparisons"][wildcards.comparison],
        pipepath = config['pipepath']
    log:
        "logs/diffexp/{comparison}.log"
    script:
        "../../scripts/diff.R"
