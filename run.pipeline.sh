#!/bin/bash


# Initialize new project
#mkdir my_project && cd my_project
#cp /scratch/group/lilab/Phil/bulkRNASeq_pipeline/run.pipeline.sh .
#cp /path/to/pipeline/template_config.yaml config.yaml
# Edit config.yaml

pipeline_dir=/scratch/group/lilab/Phil/bulkRNASeq_pipeline

# Check config
#python ${pipeline_dir}/config_validate.py config.yaml || exit 1

# Submit to SLURM
snakemake --snakefile ${pipeline_dir}/Snakefile \
    --configfile config.yaml \
    --executor slurm \
    -j 30 --latency-wait 60

