
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
        "05.Normalization_DESeq2/counts_normalized.tsv"
    output:
        "05.Normalization_DESeq2/corr.heatmap.All_vs_All.png",
        expand("05.Normalization_DESeq2/{corr}.png", sample=config['corr'])
    params:
        pipepath = config['pipepath']
    shell:
        """
python {params.pipepath}/scripts/pca.py {input} {",".join(output)}
        """

rule pca:
    input:
        "05.Normalization_DESeq2/counts_normalized.tsv"
    output:
        "05.Normalization_DESeq2/pca.plots.pdf",
        "05.Normalization_DESeq2/pca.plots.png"
    params:
        pipepath = config['pipepath']
    shell:
        """
module load Anaconda3/2024.02-1
python {params.pipepath}/scripts/pca.py {input} {",".join(output)}
        """
