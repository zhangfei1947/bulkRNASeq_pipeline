rule generate_venn:
    input:
        lambda wildcards: [
            f"06.Diff_Expression/{comp}.FC2.padj0.05.deseq2_results.tsv"
            for comp in config["venn"][wildcards.group].split(",")
        ]
    output:
        "07.Venn/{group}.venn.png"
    params:
        labels = lambda wildcards: [
            comp.strip() 
            for comp in config["venn"][wildcards.group].split(",")
        ],
        r_libs = "/scratch/group/lilab/software/R_library/4.3"
    envmodules:
        "R/4.3.1",
        "OpenMPI/4.1.4",
        "GCC/12.2.0"
    script:
        "../../scripts/venn.R"

