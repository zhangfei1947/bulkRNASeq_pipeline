#!/usr/bin/env Rscript

library(DESeq2)
library(EnhancedVolcano)

args <- commandArgs(trailingOnly = TRUE)
rc_file <- args[1]
sample <- args[2]
group <- args[3]
cmp_info <- args[4]
outpath <- args[5]

sample <- strsplit(sample, ",")[[1]]
group <- strsplit(group, ",")[[1]]
sample_info <- cbind(sample, group)
count_matrix <- read.table(rc_file, sep="\t", header=TRUE, row.names=1)

#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_info, design=~group)
#filter low readcount
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#differential expression analysis
dds <- DESeq(dds)

cmp_info = strsplit(cmp_info, ",")[[1]]

for (i in 1:length(cmp_info)){

	res <- results(dds, contrast=c("group", strsplit(cmp_info[[i]], "_vs_")[[1]]))
	resOrdered <- res[order(res$padj),]
	write.table(as.data.frame(resOrdered), file=paste0(outpath, cmp_info[[i]] ,".deseq2_results.tsv"), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

	sigs = abs(resOrdered$log2FoldChange)>1 & resOrdered$padj<0.05
	sigs[is.na(sigs)] <- FALSE
	res_sig <-  resOrdered[sigs, ]
	write.table(as.data.frame(res_sig), file=paste0(outpath, cmp_info[[i]] ,".FC2.padj0.05.deseq2_results.tsv"), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

	vocplot <- EnhancedVolcano(
		res,
	    lab=rownames(res),
	    title="", subtitle="", caption="",
	    FCcutoff=2,
	    pCutoff=0.05,
	    x='log2FoldChange', y='padj',
	    legendPosition = "bottom",
	    border = 'full',
	    )
	pdf(paste0(outpath, cmp_info[[i]] ,".volcano.pdf"))
	print(vocplot)
	dev.off()
	png(paste0(outpath, cmp_info[[i]] ,".volcano.png"))
	print(vocplot)
	dev.off()
}

