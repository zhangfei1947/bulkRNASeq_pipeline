
qc_table <- snakemake@input[[1]]
outplot <- snakemake@output[[1]]

data <- read.table(qc_table, header=TRUE, sep="\t")

data$duplication_rate <- data$duplication_rate * 100

png(outplot, width=500, height=600, res=120)

boxplot(data$duplication_rate, 
        main="Distribution of Duplication Rate",
        ylab="Duplication Rate (%)",
        ylim=c(min(data$duplication_rate)-1, max(data$duplication_rate)+2))

points(rep(1, nrow(data)), data$duplication_rate, pch=16)

highest_indices <- order(data$duplication_rate, decreasing=TRUE)[1:3]

text(rep(1.2, 3), 
     data$duplication_rate[highest_indices], 
     data$Sample[highest_indices], 
     pos=4, 
     cex=0.8)

dev.off()
