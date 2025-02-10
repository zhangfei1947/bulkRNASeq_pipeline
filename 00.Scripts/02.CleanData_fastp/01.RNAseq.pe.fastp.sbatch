#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=fastp             #Set the job name 
#SBATCH --time=1:00:00               #Set the wall clock limit
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=2          #Request 8 tasks/cores per node
#SBATCH --mem=3G                     #Request 8GB per node 
#SBATCH --output=log/fastp.%j        #Send stdout/err to "log/fastp.[jobID]"

#First Executable Line
module load GCC/11.2.0
module load fastp/0.23.2


fastp \
  --in1 $fq1 \
  --in2 $fq2 \
  --out1 $cleanfq1 \
  --out2 $cleanfq2 \
  --html ${spname}.fastp.html \
  --report_title ${spname} \
  --thread 2

echo $(date) done.main.process