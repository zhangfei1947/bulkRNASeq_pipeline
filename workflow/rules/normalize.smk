
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


def get_samples_for_groups(wildcards):
    target_groups = wildcards.groups.split("-")
    sample_list = []
    for sample_name, sample_info in config["samples"].items():
        if sample_info["group"] in target_groups:
            sample_list.append(sample_name)
    return ",".join(sample_list)

rule corr_heat:
    input:
        "05.Normalization_DESeq2/counts_normalized.tsv"
    output:
        "05.Normalization_DESeq2/corr.heatmap.All_vs_All.png",
        expand("05.Normalization_DESeq2/corr.heatmap.{groups}.png", groups=config['corr'])
    params:
        groups = config['corr']
        sps = get_samples_for_groups
        pipepath = config['pipepath']
        outpath = "05.Normalization_DESeq2"
    shell:
        """
module load Anaconda3/2024.02-1
python {params.pipepath}/scripts/corr.py {input} {";".join(params.groups)} {";".join(params.sps)} {params.outpath}
        """

rule pca:
    input:
        "05.Normalization_DESeq2/counts_normalized.tsv"
    output:
        "05.Normalization_DESeq2/pca.plots.pdf",
        "05.Normalization_DESeq2/pca.plots.png"
    params:
        color_name = config['pca_color'].keys()
        color_grp = [grp.keys() for grp in config['pca_color'].values()]
        color_sp = [grp.values() for grp in config['pca_color'].values()]
        pipepath = config['pipepath']
        outpath = "05.Normalization_DESeq2"
    shell:
        """
echo {params.color_grp}
echo {params.color_sp}
#module load Anaconda3/2024.02-1
#python {params.pipepath}/scripts/pca.py {input} 
        """
