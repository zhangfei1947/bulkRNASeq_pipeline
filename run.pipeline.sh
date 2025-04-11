#!/bin/bash

module load Anaconda3/2024.02-1
source /scratch/group/lilab/software/Python_lib/snk9_1_10/bin/activate

#export R_LIBS_USER="/scratch/group/lilab/software/R_library/4.4"
#export PATH="/scratch/group/lilab/software/jq:$PATH"
#export PATH="/scratch/group/lilab/software/subread-2.1.0/bin:$PATH"
#export PATH="/scratch/group/lilab/software/hisat2-2.2.1:$PATH"

pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline
snakemake \
 --snakefile ${pipeline_dir}/workflow/Snakefile \
 --configfile config.yaml \
 --profile ${pipeline_dir}/config