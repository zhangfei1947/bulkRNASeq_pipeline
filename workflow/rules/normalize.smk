
rule normalize_counts:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
    output:
        normalized = "05.Normalization_DESeq2/counts_normalized.tsv",
        #fpkm = "05.Normalization_DESeq2/fpkm_matrix.tsv"
    params:
#        anno = config['genome']['annotation'],
#        sample = lambda wildcards: config['samples'][wildcards.sample],
#        group = lambda wildcards: config['samples'][wildcards.sample]['group']
        r1 = lambda wildcards: config['samples'][wildcards.sample]['read1'],
        r2 = lambda wildcards: config['samples'][wildcards.sample]['read2']
    shell:
        """
echo {params.r1} >> {output.normalized}
echo {params.r2} >> {output.normalized}

#module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
#export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"
#Rscript {params.pipepath}/scripts/norm.R {input} {output}
        """