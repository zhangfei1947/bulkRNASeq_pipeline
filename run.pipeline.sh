#!/bin/bash

#snakemake module
#module load GCC/12.2.0  OpenMPI/4.1.4 snakemake/7.32.3 #deprecate
module load Anaconda3/2024.02-1
source activate snakemake732
module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.3"

pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline

# Check config
#python ${pipeline_dir}/config_validate.py config.yaml || exit 1

snakemake --snakefile ${pipeline_dir}/Snakefile \
    --configfile config.yaml \
    --profile ${pipeline_dir}/profiles/slurm
