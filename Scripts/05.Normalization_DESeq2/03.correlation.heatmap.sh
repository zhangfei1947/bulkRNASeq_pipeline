#!/bin/bash

module load Anaconda3/2024.02-1

norm_rc_file=/scratch/group/lilab/Phil/20250107_isogrp/05.Normalization_DESeq2/normalized.rc.csv

python 03.correlation.heatmap.py $norm_rc_file