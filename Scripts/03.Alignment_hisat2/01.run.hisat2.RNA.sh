#!/bin/bash

HISAT2IDX=/scratch/group/lilab/Genome_Reference/dmel_r6.55/dmel_r6.55
splicesites=/scratch/group/lilab/Genome_Reference/dmel_r6.55/dmel-all-r6.55.hisat2.splice_sites

pj_dir=YOUR_PROJECT_DIRECTORY
rawdata_dir=${pj_dir}/01.RawData
cleandata_dir=${pj_dir}/02.CleanData_fastp
mapping_dir=${pj_dir}/03.Alignment_hisat2

for spname in ${rawdata_dir}/*; 
do
	spname=$(basename "$spname")
	fq1=${cleandata_dir}/${spname}_1.fq.gz
	fq2=${cleandata_dir}/${spname}_2.fq.gz
	outdir=${mapping_dir}/${spname}
	mkdir -p $outdir
	echo $spname
	sbatch --export=outdir=$outdir,index=$HISAT2IDX,spname=$spname,fastq1=$fq1,fastq2=$fq2,splicesites=$splicesites 01.hisat2Mapping.RNA.sbatch
done
