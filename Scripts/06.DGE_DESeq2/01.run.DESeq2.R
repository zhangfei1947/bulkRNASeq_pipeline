library(DESeq2)
library(EnhancedVolcano)

args <- commandArgs(trailingOnly = TRUE)
filtered_rc_file <- args[1]
sample_info_file <- args[2]
out_path <- args[3]
setwd(out_path)

count_matrix <- read.table(filtered_rc_file,header=TRUE,row.names=1)
sample_info <- read.table(sample_info_file,header=TRUE, row.names=1)

#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_info, design=~condition)
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

mains <- c(
	"body_d1_iso/grp", 
	"body_grp_d7/d1", 
	"body_iso_d7/d1", 
	"body_d7_iso/grp", 
	"head_d1_iso/grp", 
	"head_grp_d7/d1", 
	"head_iso_d7/d1", 
	"head_d7_iso/grp"
	)

for (i in 1:length(compares)){
	res <- results(dds, contrast=strsplit(compares[[i]], "\\+")[[1]])
	resOrdered <- res[order(res$padj),]
	write.table(as.data.frame(resOrdered), file=paste0(gsub("/", ".", mains[i]),".csv"), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
	sigs = abs(resOrdered$log2FoldChange)>1 & resOrdered$padj<0.05
	sigs[is.na(sigs)] <- FALSE
	res_sig <-  resOrdered[sigs, ]
	write.table(as.data.frame(res_sig), file=paste0(gsub("/", ".", mains[i]),".logFC1.padj0.05.csv"), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
	vocplot <- EnhancedVolcano(
		res,
	    lab=rownames(res),
	    title=mains[i], subtitle="", caption="",
	    FCcutoff=2,
	    pCutoff=0.05,
	    x='log2FoldChange', y='padj',
	    legendPosition = "bottom",
	    border = 'full',
	    )
	pdf(paste0(gsub("/", ".", mains[i]),".pdf"))
	print(vocplot)
	dev.off()
	png(paste0(gsub("/", ".", mains[i]),".png"))
	print(vocplot)
	dev.off()
}




