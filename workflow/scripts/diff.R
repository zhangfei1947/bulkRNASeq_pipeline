#!/usr/bin/env Rscript

library(DESeq2)
library(EnhancedVolcano)

rc_file <- snakemake@input$counts
sample <- snakemake@params$sample
group <- snakemake@params$group
cmp_info <- snakemake@params$comparison
gene_info <- snakemake@params$anno
outpath <- snakemake@params$outpath

sample <- strsplit(sample, ",")[[1]]
group <- strsplit(group, ",")[[1]]
sample_info <- cbind(sample, group)
gene_info <- read.table(file=gene_info, header=TRUE, sep=",", row.names=1)
#create gene symbol mapping
mapping <- gene_info[, 1, drop = TRUE]
names(mapping) <- rownames(gene_info)

count_matrix <- read.table(rc_file, sep="\t", header=TRUE, row.names=1)
count_matrix <- round(count_matrix)

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
	# mapping gene symbol
	geneIDs <- rownames(res)
	genesymbols <- mapping[geneIDs]
	print(genesymbols)
	genesymbols <- ifelse(is.na(genesymbols), geneIDs, genesymbols)

	vocplot <- EnhancedVolcano(
		res,
	    lab=genesymbols,
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

