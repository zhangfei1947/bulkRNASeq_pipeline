#!/bin/bash


# Initialize new project
#mkdir my_project && cd my_project
#cp /path/to/pipeline/template_config.yaml config.yaml
# Edit config.yaml

module load GCC/12.3.0  OpenMPI/4.1.5 snakemake/8.4.2

pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline

# Check config
#python ${pipeline_dir}/config_validate.py config.yaml || exit 1

# Submit to SLURM
snakemake --snakefile ${pipeline_dir}/Snakefile \
    --configfile config.yaml \
    --cluster "sbatch -t {resources.runtime} --mem {resources.mem_mb} -c {threads} -J {rule} -o logs/slurm/{rule}_{wildcards} -e logs/slurm/{rule}_{wildcards}" \
    -j 30 --latency-wait 60

