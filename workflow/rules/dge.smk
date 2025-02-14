rule differential_expression:
    input:
        "04.Quant_featureCounts/counts_raw.tsv"
    output:
        "06.DGE_DESeq2/{comparison}_results.tsv"
    params:
        contrast = lambda wildcards: config['diff_comparisons'][wildcards.comparison],
        pipepath = config['pipepath']
    log:
        "logs/diff/{comparison}.log"
    shell:
        """
module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"
Rscript {params.pipepath}/scripts/diff.R {input} {params.contrast}
        """