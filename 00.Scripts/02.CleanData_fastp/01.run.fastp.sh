#!/bin/bash

rawdata_dir=/scratch/group/lilab/Phil/20250107_isogrp/01.RawData
cleandata_dir=/scratch/group/lilab/Phil/20250107_isogrp/02.CleanData
for spname in ${rawdata_dir}/*_*_*; do
	spname=$(basename "$spname") 
	fq1=$(ls -d ${rawdata_dir}/${spname}/*_1.fq.gz 2>/dev/null)
	fq2=$(ls -d ${rawdata_dir}/${spname}/*_2.fq.gz 2>/dev/null)
	cleanfq1=${cleandata_dir}/${spname}_1.fq.gz
	cleanfq2=${cleandata_dir}/${spname}_2.fq.gz
	echo $spname
	sbatch --export=fq1=$fq1,fq2=$fq2,cleanfq1=$cleanfq1,cleanfq2=$cleanfq2,spname=$spname 01.RNAseq.pe.fastp.sbatch
done
