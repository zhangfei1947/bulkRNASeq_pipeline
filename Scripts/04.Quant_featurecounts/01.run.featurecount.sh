#!/bin/bash

refgtf=/scratch/group/lilab/Genome_Reference/dmel_r6.55/dmel-all-r6.55.gtf
mapping_dir=/scratch/group/lilab/Phil/20250107_isogrp/03.Alignment
quant_dir=/scratch/group/lilab/Phil/20250107_isogrp/04.Quantification

for spname in ${mapping_dir}/*_*; 
do
    spname=$(basename "$spname")
    sortbam=${mapping_dir}/${spname}/${spname}.sorted.bam
    outdir=${quant_dir}/${spname}
    mkdir -p $outdir
    echo $spname
    sbatch --export=outdir=$outdir,refgtf=$refgtf,sortbam=$sortbam,spname=$spname 01.featurecount.sbatch
done
