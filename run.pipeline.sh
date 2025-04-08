#!/bin/bash

#snakemake module
#module load GCC/12.2.0  OpenMPI/4.1.4 snakemake/7.32.3 #deprecate

module load Anaconda3/2024.02-1
source activate snakemake913


pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline

# Check config
#python ${pipeline_dir}/config_validate.py config.yaml || exit 1

snakemake \
  --snakefile ${pipeline_dir}/workflow/Snakefile \
  --configfile config.yaml \
  --executor slurm \
  --use-singularity

