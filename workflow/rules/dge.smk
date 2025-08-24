

rule diff_analysis:
    input:
        counts = "04.Quant_featureCounts/counts_filter.tsv"
    output:
        results = expand("06.Diff_Expression/{comparison}.deseq2_results.tsv", comparison=config['diff_comparisons']),
        sigresults = expand("06.Diff_Expression/{comparison}.FC2.padj0.05.deseq2_results.tsv", comparison=config['diff_comparisons'])
    params:
        comparison = ",".join(config["diff_comparisons"]),
        sample = ",".join( config['samples'].keys() ),
        group = ",".join( sample['group'] for sample in config["samples"].values() ),
        anno = config['genome']['geneloc'],
        outpath = "06.Diff_Expression/",
    script:
        "../scripts/diff.R" 
