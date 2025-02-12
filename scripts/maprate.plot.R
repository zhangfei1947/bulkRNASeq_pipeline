library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
map_table <- args[1]
outplot <- args[2]

data <- read.table(map_table, header=TRUE, sep = "\t")
# Create boxplot
p <- ggplot(data, aes(x="", y=Exact.Match.Rate)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "Boxplot of Exact Match Rate",
       x = "",
       y = "Exact Match Rate (%)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
# Add text labels for the three lowest points
outliers <- subset(data, Exact.Match.Rate %in% c(min(Exact.Match.Rate, na.rm = TRUE), 
                                               min(Exact.Match.Rate, na.rm = TRUE),
                                               min(Exact.Match.Rate, na.rm = TRUE)))
p <- p + geom_text(aes(label=Sample.Name, y=Exact.Match.Rate), vjust=-0.5)  
# Save the plot as PNG
ggsave(outplot, plot=p, width=6, height=6, dpi=300)
