#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=hisat2             #Set the job name 
#SBATCH --time=2:00:00               #Set the wall clock limit
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=4          #Request 8 tasks/cores per node
#SBATCH --mem=12G                     #Request xGB per node 
#SBATCH --output=log/%x.%j        #Send stdout/err to "log/xxx.[jobID]"

#First Executable Line
module load Anaconda3/2024.02-1
source activate hisat2

##output bam sorted by read name, for the featurecounts

hisat2 \
 -x $index \
 -1 $fastq1 -2 $fastq2 \
 --rna-strandness FR \
 --known-splicesite-infile $splicesites \
 --summary-file ${outdir}/${spname}.summary \
 --no-unal \
 -p 4 \
 | samtools view -bS 2> ${outdir}/samview.log | samtools sort -o ${outdir}/${spname}.sorted.bam 2> ${outdir}/samsort.log

 echo $(date) done.main.process