#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=rna_snake             #Set the job name 
#SBATCH --time=4:00:00               #Set the wall clock limit
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=32          #Request 8 tasks/cores per node
#SBATCH --mem=60G                     #Request xGB per node 
#SBATCH --output=log/%x.%j        #Send stdout/err to "log/xxx.[jobID]"

#First Executable Line
module load Anaconda3/2024.02-1
source activate snakemake

pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline

# Check config
#python ${pipeline_dir}/config_validate.py config.yaml || exit 1

snakemake --snakefile ${pipeline_dir}/Snakefile \
    --configfile config.yaml \
    -j 15 --latency-wait 60
