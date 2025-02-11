#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=snk             #Set the job name
#SBATCH --time=4:00:00               #Set the wall clock limit
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=1          #Request 8 tasks/cores per node
#SBATCH --mem=8G                     #Request 8GB per node
#SBATCH --output=snk.%j        #Send stdout/err to "log/fastp.[jobID]"


#snakemake module
module load GCC/12.2.0  OpenMPI/4.1.4 snakemake/7.32.3

pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline

# Check config
#python ${pipeline_dir}/config_validate.py config.yaml || exit 1

snakemake --snakefile ${pipeline_dir}/Snakefile \
    --configfile config.yaml \
    --cluster "sbatch -t {resources.runtime} --mem {resources.mem_mb} -c {threads} -N {resources.nodes} -n {resources.ntasks_per_node} -J {rule}" \
    -j 30 --latency-wait 60 --rerun-incomplete
