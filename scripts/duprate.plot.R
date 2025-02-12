library(ggplot2)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
qc_table <- args[1]


data <- read.delim(qc_table, header = TRUE, sep = "\t")


data <- data %>%
  mutate(rank = rank(-duplication_rate)) %>%
  filter(rank <= 3)


ggplot(data, aes(x = "", y = duplication_rate, fill = Sample)) +
  geom_boxplot(width = 0.5) +
  geom_text(aes(label = Sample), vjust = -1, position = position_jitter(width = 0.2, height = 0)) +
  theme_minimal() +
  labs(title = "Boxplot of Duplication Rate",
       x = "",
       y = "Duplication Rate") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("duprate.boxplot.png", width=6, height=6, dpi=300)
