#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
filtered_rc_file <- args[1]
sample_info_file <- args[2]
out_path <- args[3]
anno_file <- args[4]

input_file <- snakemake@input[[1]]
output_norm <- snakemake@output$normalized
output_fpkm <- snakemake@output$fpkm

.libPaths(c(snakemake@params$r_libs, .libPaths()))

library(DESeq2)

count_matrix <- read.table(filtered_rc_file,header=TRUE,row.names=1)
sample_info <- read.table(sample_info_file,header=TRUE, row.names=1)

#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_info, design=~condition)

#median of ratios method of normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file=paste0(out_path,"/counts_normalized.tsv"), sep="\t", quote=F, col.names=T, row.names=T)

#cal FPKM
library(edgeR)

gene_info = read.csv(anno_file, header=T, row.names=1)
gene_info = gene_info[row.names(normalized_counts),]

genefpkm <- rpkm(normalized_counts, gene.length=gene_info$Gene_length, log=FALSE)
write.table(genefpkm, file = paste0(out_path,"/fpkm_matrix.tsv"), sep="\t", quote=F, row.names=T)