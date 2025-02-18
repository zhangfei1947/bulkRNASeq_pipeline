localrules: corr_heat

rule normalize_counts:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
    output:
        normalized = "05.Normalization_DESeq2/counts_normalized.tsv",
        fpkm = "05.Normalization_DESeq2/fpkm_matrix.tsv"
    params:
        anno = config['genome']['geneloc'],
        sample = ",".join( config['samples'].keys() ),
        group = ",".join( sample['group'] for sample in config["samples"].values() ),
        pipepath = config['pipepath']
    shell:
        """
module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"
Rscript {params.pipepath}/scripts/norm.R {input} {output.normalized} {output.fpkm} {params.anno} {params.sample} {params.group}
        """


rule corr_heat:
    input:
        norm_counts = "05.Normalization_DESeq2/counts_normalized.tsv"
    output:
        "05.Normalization_DESeq2/corr.heatmap.{corr_name}.png"
    params:
        target_groups = lambda wildcards: config["corr"][wildcards.corr_name].split(","),
        sample_mapping = lambda wildcards: (print(f"Current corr_name: {wildcards.corr_name}"), {s: info["group"] for s, info in config["samples"].items()})[-1]
        pipepath = config['pipepath']
    script:
        """
{params.pipepath}/scripts/corr.py
        """

#rule pca:
#    input:
#        "05.Normalization_DESeq2/counts_normalized.tsv"
#    output:
#        "05.Normalization_DESeq2/pca.plots.pdf",
#        "05.Normalization_DESeq2/pca.plots.png"
#    params:
#        color_name = config['pca_color'].keys(),
#        color_grp = [grp.keys() for grp in config['pca_color'].values()],
#        color_sp = [grp.values() for grp in config['pca_color'].values()],
#        pipepath = config['pipepath'],
#        outpath = "05.Normalization_DESeq2"
#    script:
#        """
#{params.pipepath}/scripts/pca.py
#        """
#