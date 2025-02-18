#!/usr/bin/env Rscript

library(DESeq2)
library(EnhancedVolcano)

counts_file <- snakemake@input$counts
output_files <- snakemake@output
sample_mapping <- snakemake@params$sample_info
cmp_info <- snakemake@params$comparison
log_file <- snakemake@log[[1]]
sink(log_file)

print(output_files)
print(sample_mapping)
print(cmp_info)


count_matrix <- read.table(counts_file, sep="\t", header=TRUE, row.names=1)

#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_info, design=~group)
#filter low readcount
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#differential expression analysis
dds <- DESeq(dds)

compares <- c(
	"condition+d1_iso_b+d1_grp_b",
	"condition+d7_grp_b+d1_grp_b",
	"condition+d7_iso_b+d1_iso_b",
	"condition+d7_iso_b+d7_grp_b",

	"condition+d1_iso_h+d1_grp_h",
	"condition+d7_grp_h+d1_grp_h",
	"condition+d7_iso_h+d1_iso_h",
	"condition+d7_iso_h+d7_grp_h"
	)

for (i in 1:length(compares)){
	res <- results(dds, contrast=strsplit(compares[[i]], "\\+")[[1]])
	resOrdered <- res[order(res$padj),]
	write.table(as.data.frame(resOrdered), file=output_files$results, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

	sigs = abs(resOrdered$log2FoldChange)>1 & resOrdered$padj<0.05
	sigs[is.na(sigs)] <- FALSE
	res_sig <-  resOrdered[sigs, ]
	write.table(as.data.frame(res_sig), file=output_files$sigresults, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

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
	pdf(output_files$vcpdf)
	print(vocplot)
	dev.off()
	png(output_files$vcpng)
	print(vocplot)
	dev.off()
}

