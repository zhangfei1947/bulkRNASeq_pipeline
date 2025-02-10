#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=featureCounts     #Set the job name 
#SBATCH --time=00:30:00               #Set the wall clock limit
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=4          #Request 8 tasks/cores per node
#SBATCH --mem=8G                    #Request xGB per node 
#SBATCH --output=log/%x.%j           #Send stdout/err to "log/xxx.[jobID]"

#First Executable Line
module load GCC/12.3.0 Subread/2.0.8

featureCounts \
  -a $refgtf \
  -o ${outdir}/${spname}.fC \
  -F GTF \
  -t exon \
  -g gene_id \
  -O -M --fraction \
  -s 2 \
  -p \
  --countReadPairs \
  -P \
  -B \
  -d 40 \
  -D 800 \
  -T 4 \
  $sortbam

 echo $(date) done.main.process
