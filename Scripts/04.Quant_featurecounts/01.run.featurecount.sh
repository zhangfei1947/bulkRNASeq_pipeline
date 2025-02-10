#!/bin/bash

refgtf=/scratch/group/lilab/Genome_Reference/dmel_r6.55/dmel-all-r6.55.gtf

pj_dir=YOUR_PROJECT_DIRECTORY
rawdata_dir=${pj_dir}/01.RawData
mapping_dir=${pj_dir}/03.Alignment_hisat2
quant_dir=${pj_dir}/04.Quant_featurecounts

for spname in ${rawdata_dir}/*; 
do
    spname=$(basename "$spname")
    sortbam=${mapping_dir}/${spname}/${spname}.sorted.bam
    outdir=${quant_dir}/${spname}
    mkdir -p $outdir
    echo $spname
    sbatch --export=outdir=$outdir,refgtf=$refgtf,sortbam=$sortbam,spname=$spname 01.featurecount.sbatch
done
