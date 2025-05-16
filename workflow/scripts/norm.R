#!/usr/bin/env Rscript

rc_file <- snakemake@input[[1]]
anno_file <- snakemake@params$anno
sample <- snakemake@params$sample
group <- snakemake@params$group
output_norm <- snakemake@output$normalized
output_fpkm <- snakemake@output$fpkm
output_vst <- snakemake@output$vst


library(DESeq2)

count_matrix <- read.table(rc_file, header=TRUE, row.names=1)
count_matrix <- round(count_matrix)
sample <- strsplit(sample, ",")[[1]]
groups <- strsplit(group, ",")[[1]]
sample_info <- cbind(sample, groups)

print(sample_info)
#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_info, design=~groups)

#median of ratios method of normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file=paste0(output_norm), sep="\t", quote=F, col.names=T, row.names=T)

#vst normalization for heatmap plot
vsd <- vst(dds)
write.table(round(assay(vsd),2), file=output_vst, sep="\t", quote=F, col.names=T, row.names=T)

#cal FPKM
library(edgeR)

gene_info = read.csv(anno_file, header=T, row.names=1)
gene_info = gene_info[row.names(normalized_counts),]

genefpkm <- rpkm(normalized_counts, gene.length=gene_info$Gene_length, log=FALSE)
write.table(genefpkm, file = paste0(output_fpkm), sep="\t", quote=F, row.names=T)