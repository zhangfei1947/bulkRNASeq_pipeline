#!/bin/bash

module restore R431
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"

filtered_rc_file="/scratch/group/lilab/Phil/20250107_isogrp/04.Quantification/merged_featurecounts_filtered.tsv"
out_path="/scratch/group/lilab/Phil/20250107_isogrp/06.DGE_DESeq2"
sample_info_file="/scratch/group/lilab/Phil/20250107_isogrp/05.Normalization_DESeq2/sampleinfo.csv"

Rscript 01.run.DESeq2.R $filtered_rc_file $sample_info_file $out_path 