#!/bin/bash

module load Anaconda3/2024.02-1
source activate snakemake732

export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"

pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline


snakemake --snakefile ${pipeline_dir}/Snakefile --configfile config.yaml 

