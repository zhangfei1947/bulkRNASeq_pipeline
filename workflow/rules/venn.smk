rule generate_venn:
    input:
        lambda wildcards: [
            f"06.Diff_Expression/{comp}.FC2.padj0.05.deseq2_results.tsv"
            for comp in config["venn"][wildcards.group].split(",")
        ]
    output:
        "07.Venn/{group}.venn.png"
    params:
        input_files = lambda wildcards: ",".join([
            f"06.Diff_Expression/{comp.strip()}.FC2.padj0.05.deseq2_results.tsv"
            for comp in config["venn"][wildcards.group].split(",")
            if comp.strip()
        ]),
        labels_str = lambda wildcards: ",".join([
            comp.strip()
            for comp in config["venn"][wildcards.group].split(",")
            if comp.strip()
        ])
    script:
        "../scripts/venn.R"
