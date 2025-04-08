#!/bin/bash

module load Anaconda3/2024.02-1
source activate snakemake913

pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline

snakemake \
  --snakefile ${pipeline_dir}/workflow/Snakefile \
  --configfile config.yaml \
  --executor slurm \
  --use-singularity

