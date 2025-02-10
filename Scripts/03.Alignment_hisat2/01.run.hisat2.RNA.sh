#!/bin/bash

HISAT2IDX=/scratch/group/lilab/Genome_Reference/dmel_r6.55/dmel_r6.55
splicesites=/scratch/group/lilab/Genome_Reference/dmel_r6.55/dmel-all-r6.55.hisat2.splice_sites
rawdata_dir=/scratch/group/lilab/Phil/20250107_isogrp/01.RawData
cleandata_dir=/scratch/group/lilab/Phil/20250107_isogrp/02.CleanData
mapping_dir=/scratch/group/lilab/Phil/20250107_isogrp/03.Alignment


for spname in ${rawdata_dir}/*_*; 
do
	spname=$(basename "$spname")
	fq1=${cleandata_dir}/${spname}_1.fq.gz
	fq2=${cleandata_dir}/${spname}_2.fq.gz
	outdir=${mapping_dir}/${spname}
	mkdir -p $outdir
	echo $spname
	sbatch --export=outdir=$outdir,index=$HISAT2IDX,spname=$spname,fastq1=$fq1,fastq2=$fq2,splicesites=$splicesites 01.hisat2Mapping.RNA.sbatch
done
