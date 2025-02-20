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
        pipepath = config['pipepath']
    shell:
        """
module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"
Rscript {params.pipepath}/scripts/venn.R {input} {params.labels} {output}
        """
