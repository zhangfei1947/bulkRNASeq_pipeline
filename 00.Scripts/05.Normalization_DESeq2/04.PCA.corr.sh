#!/bin/bash

module load Anaconda3/2024.02-1

norm_rc_file=/scratch/group/lilab/Phil/20250107_isogrp/05.Normalization_DESeq2/normalized.rc.csv
outpath=/scratch/group/lilab/Phil/20250107_isogrp/05.Normalization_DESeq2

python 04.PCA.py $norm_rc_file $outpath

python 04.sample.corr.py $norm_rc_file $outpath
