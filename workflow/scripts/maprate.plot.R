
args <- commandArgs(trailingOnly = TRUE)
map_table <- snakemake@input[[1]]
outplot <- snakemake@output[[1]]

png(outplot, width=500, height=600, res=120)

data <- read.table(map_table, header=TRUE, sep="\t")

boxplot(data$Exact.Match.Rate, 
        main="Distribution of Exact Match Rate",
        ylab="Exact Match Rate (%)",
        ylim=c(min(data$Exact.Match.Rate)-0.5, max(data$Exact.Match.Rate)+0.5))

points(rep(1, nrow(data)), data$Exact.Match.Rate, pch=16)

lowest_indices <- order(data$Exact.Match.Rate)[1:3]

text(rep(1.2, 3), 
     data$Exact.Match.Rate[lowest_indices], 
     data$Sample.Name[lowest_indices], 
     pos=4, 
     cex=0.8)

dev.off()
