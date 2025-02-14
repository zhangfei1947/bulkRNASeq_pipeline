#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

rc_file <- args[1]
anno_file <- args[4]
sample <- args[5]
group <- args[6]
output_norm <- args[2]
output_fpkm <- args[3]

library(DESeq2)

count_matrix <- read.table(rc_file, skip=1, header=TRUE, row.names=1)
samples <- strsplit(sample, ",")[[1]]
groups <- strsplit(group, ",")[[1]]
sample_info <- cbind(samples, groups)

#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_info, design=~groups)

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