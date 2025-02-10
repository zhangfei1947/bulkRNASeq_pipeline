library(VennDiagram)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
sig_files <- args[1]
cat_names <- args[2]
out_path <- args[3]
setwd(out_path)


draw_venn <- function(files, category_names, output_prefix="venn_diagram") {
  
  # Validate input length
  if (length(files) < 2 || length(files) > 4) {
    stop("This function supports 2 to 4 input files.")
  }
  if (length(files) != length(category_names)) {
    stop("Number of files must match number of category names.")
  }
  
  # Read gene lists from files
  gene_lists <- list()
  for (i in seq_along(files)) {
    # Read first column (gene names) from each file
    data <- read.table(files[i], header=TRUE, sep="\t")
    gene_lists[[category_names[i]]] <- as.character(row.names(data))
  }
  
  # Set Venn diagram parameters
  venn_params <- list(
    x = gene_lists,
    category.names = category_names,
    col = rainbow(length(gene_lists)),  # Colorful outline
    fill = NA,                         # No fill color
    alpha = 0.3,                       # Transparency (not used with fill=NA)
    lwd = 2,                           # Line width
    cex = 1.2,                         # Label size
    cat.cex = 1.2,                     # Category name size
    cat.pos = c(0, 0, 0)[1:length(gene_lists)],  # Category position
    cat.dist = c(0.05, 0.05, 0.05)[1:length(gene_lists)]  # Category distance
  )
  
  # Generate output files
  output_files <- c(
    svg = paste0(output_prefix, ".svg"),
    png = paste0(output_prefix, ".png")
  )
  
  # Create svg
  do.call(venn.diagram, c(venn_params, list(
    filename = output_files["svg"],
    imagetype = "svg"
  )))
  
  # Create PNG (300 dpi)
  do.call(venn.diagram, c(venn_params, list(
    filename = output_files["png"],
    imagetype = "png",
    height = 1200,
    width = 1200,
    resolution = 300
  )))
  
  message("Venn diagrams saved to:\n", paste(output_files, collapse="\n"))
}


files <- strsplit(sig_files, ",")[[1]]  
category_names <- strsplit(cat_names, ",")[[1]]  

draw_venn(files, category_names, paste(category_names, collapse="_"))
