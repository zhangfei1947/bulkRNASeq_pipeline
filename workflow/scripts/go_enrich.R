library(clusterProfiler)
library(org.Dm.eg.db)
library(ggplot2)
 
input_file <- snakemake@input$diff
output_file <- snakemake@output$res
output_pdf <- snakemake@output$pdf
output_png <- snakemake@output$png

gene_list <- row.names(read.table(input_file, header=TRUE, sep="\t"))

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
  if (ont=="ALL" && !is.null(result) && length(result) > 0) {
    write.csv(result, output_file, row.names=FALSE)
    plot <- dotplot(result, showCategory=20, font.size=9) +
      theme(plot.title=element_blank(),
            axis.title.x=element_blank(),
            plot.margin=margin(5.5, 40, 5.5, 5.5, "pt")) +
      annotate("text", x=Inf, y=2, label=ont, hjust=1.5, vjust=1.5, size=6)
    ggsave(filename=output_pdf, plot=plot, width=8, height=9)
    ggsave(filename=output_png, plot=plot, width=8, height=9, dpi=300)
  }
}