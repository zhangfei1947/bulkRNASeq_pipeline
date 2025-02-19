
args <- commandArgs(trailingOnly = TRUE)
fc_table <- args[1]
outplot <- args[2]

data <- read.table(fc_table, header=TRUE, sep="\t")

colnames(data)[-1] <- sapply(strsplit(colnames(data)[-1], "[.]"), function(x) x[3])

data <- data[rowSums(data[,-1]) > 0,]

plot_data <- as.matrix(t(data[,-1]))

png_width <- ncol(plot_data) * 50 + 500

png(outplot, width=png_width, height=600, res=120)

colors <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C")

# Create stacked barplot

barplot(t(plot_data),
        col=colors[1:nrow(data)],
        main="Read Assignment Distribution",
        legend.text=data$Status,
        args.legend=list(x="bottom",inset=c(0,0.1), cex=0.7, bg="white", box.lty=0),
        ylab="Number of Reads",
        las=2)  

dev.off()
