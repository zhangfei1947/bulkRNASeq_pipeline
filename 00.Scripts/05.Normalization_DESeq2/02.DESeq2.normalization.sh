#!/bin/bash

module restore R431
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"

filtered_rc_file=/scratch/group/lilab/Phil/20250107_isogrp/04.Quantification/merged_featurecounts_filtered.tsv
outdir=/scratch/group/lilab/Phil/20250107_isogrp/05.Normalization_DESeq2
sample_info_file=${outdir}/sampleinfo.csv
gene_info=/scratch/group/lilab/Genome_Reference/dmel_r6.55/gene_annotation.csv


Rscript 02.DESeq2.normalization.R $filtered_rc_file $sample_info_file $outdir $gene_info
