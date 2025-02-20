rule go_enrichment:
    input:
        diff = "06.Diff_Expression/{comparison}.FC2.padj0.05.deseq2_results.tsv"
    output:
        outfile = "08.Enrichment/{comparison}.go_enrichment.tsv",
        pipepath = config['pipepath']
    shell:
        """
module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"
Rscript {params.pipepath}/scripts/go_enrich.R {input.diff} {output.outfile}
        """
