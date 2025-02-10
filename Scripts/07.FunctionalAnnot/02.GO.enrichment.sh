#!/bin/bash

module restore R431
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"

output_dir=/scratch/group/lilab/Phil/20250107_isogrp/07.FunctionalAnnot

for csv in $(find ${output_dir} -name "*logFC1.padj0.05.csv"); do
	echo $csv
	input_file=$csv
    output_prefix=${csv%.logFC1.padj0.05.csv}
    Rscript 02.GO.enrichment.R $input_file $output_prefix $output_dir
done
