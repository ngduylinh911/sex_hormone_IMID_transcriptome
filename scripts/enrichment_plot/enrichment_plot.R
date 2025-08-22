library(ggplot2)  # for plotting
library(dplyr)   # for data manipulation

rm(list = ls()) # Clear the workspace
setwd("C:/Users/ngduy/Documents/Project/Final_code/code/enrichment_plot")

#--- Plotting KEGG enrichment results for estrogen DEGs ---
# Load the KEGG enrichment results for estrogen
estrogen_kegg_enrich <- read.csv("KEGG_2021_Human_table_estrogen.txt", header = TRUE, sep = "\t")
estrogen_kegg_enrich <- estrogen_kegg_enrich[order(estrogen_kegg_enrich$Adjusted.P.value),]
estrogen_kegg_enrich$log_FDR <- -log10(estrogen_kegg_enrich$Adjusted.P.value)

# Optional: select top terms
top_estrogen_enriched <- estrogen_kegg_enrich[1:10, ]

# Create barplot
ggplot(top_estrogen_enriched, aes(x = reorder(Term, log_FDR), y = log_FDR)) +
  geom_col(fill = "lightsalmon", color = "grey") + 
  coord_flip() +
  labs(
    title = "Top Enriched KEGG Pathways for estrogen DEGs",
    x = "Pathways",
    y = "Adjusted P-value (-log10 scale)"
  ) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_text(size = 15, hjust = 0),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))  

# Save the plot
ggsave("../../result/pbmc/img/estrogen_kegg_enrichment_plot.jpg", width = 12, height = 5, dpi = 300)

#--- End of KEGG enrichment results for estrogen DEGs ---

#--- Plotting KEGG enrichment results for testosterone DEGs ---

# Load the KEGG enrichment results for testosterone
testosterone_kegg_enrich <- read.csv("KEGG_2021_Human_table_testosterone.txt", header = TRUE, sep = "\t")
testosterone_kegg_enrich <- testosterone_kegg_enrich[order(testosterone_kegg_enrich$Adjusted.P.value),]
testosterone_kegg_enrich$log_FDR <- -log10(testosterone_kegg_enrich$Adjusted.P.value)

# Optional: select top terms
top_testosterone_enriched <- testosterone_kegg_enrich[1:10, ]

# Create barplot
ggplot(top_testosterone_enriched, aes(x = reorder(Term, log_FDR), y = log_FDR)) +
  geom_col(fill = "lightsalmon", color = "grey") + 
  coord_flip() +
  labs(
    title = "Top Enriched KEGG Pathways for testosterone DEGs",
    x = "Pathways",
    y = "Adjusted P-value (-log10 scale)"
  ) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_text(size = 15, hjust = 0),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))  

# Save the plot
ggsave("../../result/pbmc/img/testosterone_kegg_enrichment_plot.jpg", width = 12, height = 5, dpi = 300)

#--- End of KEGG enrichment results for testosterone DEGs ---
