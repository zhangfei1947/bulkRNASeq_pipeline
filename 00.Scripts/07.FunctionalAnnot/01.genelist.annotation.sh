#!/bin/bash

module load Anaconda3/2024.02-1

anno_file=/scratch/group/lilab/Genome_Reference/dmel_r6.55/gene_annotations.tsv
input_dir=/scratch/group/lilab/Phil/20250107_isogrp/06.DGE_DESeq2
postfix=padj0.05.csv
output_dir=/scratch/group/lilab/Phil/20250107_isogrp/07.FunctionalAnnot


python 01.genelist.annotation.py  $anno_file  $input_dir  $postfix  $output_dir
