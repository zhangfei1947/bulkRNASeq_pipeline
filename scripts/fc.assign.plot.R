
args <- commandArgs(trailingOnly = TRUE)
fc_table <- args[1]
outplot <- args[2]

data <- read.table(fc_table, header=TRUE, sep="\t")

colnames(data)[-1] <- sapply(strsplit(colnames(data)[-1], "[.]"), function(x) x[3])

data <- data[rowSums(data[,-1]) > 0,]

plot_data <- as.matrix(t(data[,-1]))

png_width <- ncol(plot_data) * 160 + 300

png(outplot, width=png_width, height=600, res=120)

colors <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C")

# Create stacked barplot
par(xpd=TRUE)
barplot(t(plot_data),
        col=colors[1:nrow(data)],
        legend.text=data$Status,
        args.legend=list(x="right",inset=c(-1,0), cex=0.8, bg="white", box.lty=0),
        main="Read Assignment Distribution",
        ylab="Number of Reads",
        las=2)  

dev.off()
