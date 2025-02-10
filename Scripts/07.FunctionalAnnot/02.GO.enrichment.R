library(clusterProfiler)
library(org.Dm.eg.db)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_prefix <- args[2]
output_dir <- args[3]
setwd(output_dir)

gene_list <- read.table(input_file, header=TRUE, sep="\t", quote='"')[,1]

enrichment_results <- list(
  BP = enrichGO(gene          = gene_list,
                OrgDb         = org.Dm.eg.db,
                keyType       = "FLYBASE",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2),
  MF = enrichGO(gene          = gene_list,
                OrgDb         = org.Dm.eg.db,
                keyType       = "FLYBASE",
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2),
  CC = enrichGO(gene          = gene_list,
                OrgDb         = org.Dm.eg.db,
                keyType       = "FLYBASE",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2),
  ALL = enrichGO(gene          = gene_list,
                OrgDb         = org.Dm.eg.db,
                keyType       = "FLYBASE",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2)
)

# Plot and save the top 10 enriched GO terms for each ontology
for (ont in names(enrichment_results)) {
  result <- enrichment_results[[ont]]
  if (!is.null(result) && length(result) > 0) {
    write.csv(result, paste0(output_prefix, ".GO_enrich.csv"), row.names=FALSE)
    plot <- dotplot(result, showCategory = 20, font.size=9) +
      theme(plot.title = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(5.5, 40, 5.5, 5.5, "pt")) +
      annotate("text", x = Inf, y = 2, label = ont, hjust = 1.5, vjust = 1.5, size = 6)
    ggsave(filename = paste0(output_prefix, "_", ont, "_top20_enrichment.png"), plot = plot, width = 8, height = 9, dpi = 100)
  }
}