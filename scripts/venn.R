library(VennDiagram)

#input_files <- snakemake@input[[1]]    
#labels <- snakemake@params$labels     
#output_file <- snakemake@output[[1]]   

args <- commandArgs(trailingOnly=TRUE)
input_files <- strsplit(args[1], ",")[[1]]
labels <- strsplit(args[2], ",")[[1]]
output_file <- args[3]


read_diff_genes <- function(file) {
  df <- read.delim(file, sep = "\t")
#  sig_genes <- df[
#    df$padj < 0.05 & 
#    abs(df$log2FoldChange) >= 1, 
#    "gene_id"  
#  ]
  sig_genes <- row.names(df)
  return(unique(sig_genes)) 
}

gene_lists <- lapply(input_files, read_diff_genes)
names(gene_lists) <- labels  

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

venn.diagram(
  x = gene_lists,
  category.names = names(gene_lists),
  filename = output_file,
  imagetype = "png",
  height = 3000, 
  width = 3000,
  resolution = 300,
  col = "black",
  fill = c("#E69F00", "#56B4E9", "#009E73", "#F0E442")[1:length(gene_lists)],
  alpha = 0.5,
  cex = 2.5,
  fontfamily = "sans",
  cat.cex = 2.2,
  cat.fontfamily = "sans",
  margin = 0.1
)
