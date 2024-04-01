# written by Crystal Shan 04/2024
library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)

table <- data.table::fread("~/Downloads/combined_gene_expression.csv")
gene_of_interest <- table[1,]$Gene

plot_gene_expression <- function(expression_table, gene_of_interest) {
  
  data <- expression_table
  
  # Filtering the data for the gene of interest
  filtered_data <- subset(data, Gene == gene_of_interest)
  
  # Creating the box-whisker plot
  p <- ggplot(filtered_data, aes(x=Species, y=Expression, fill=Species)) +
    geom_boxplot() +
    labs(title=paste("Expression of", gene_of_interest, "by Species"),
         x="Species", y="Expression") +
    theme_minimal() +
    scale_fill_brewer(palette="Set1")
  
  print(p)
}

save.image(file="temp_test.RData")
