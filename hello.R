# Package: Fun
# Description: A comprehensive R package for GO & KEGG enrichment analysis with publication-ready visualizations.
# Uniqueness: Combines GO & KEGG in a single figure, automates the workflow, and ensures robust error handling.

# Load required packages
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)
library(ggpubr)

# Function to perform GO and KEGG analysis with ggplot visualization
perform_go_kegg_analysis <- function(file) {
  # Read input gene list
  df <- read.csv(file)
  if (!"Gene_Name" %in% colnames(df)) stop("Error: Gene_Name column missing in CSV file")

  # Convert gene symbols to Entrez IDs
  gene_list <- bitr(df$Gene_Name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(gene_list) == 0) stop("Error: No valid gene mappings found")

  # Prepare geneList for GSEA
  df2 <- merge(df, gene_list, by.x = "Gene_Name", by.y = "SYMBOL")
  gene_list_sorted <- df2$log2FoldChange
  names(gene_list_sorted) <- df2$ENTREZID
  gene_list_sorted <- sort(gene_list_sorted, decreasing = TRUE)

  # Perform GO enrichment analysis
  go_results <- gseGO(geneList = gene_list_sorted, ont = "ALL", keyType = "ENTREZID", OrgDb = org.Hs.eg.db,
                      pvalueCutoff = 0.1, verbose = TRUE, pAdjustMethod = "BH")
  if (is.null(go_results)) stop("Error: No GO enrichment results found")

  # Perform KEGG enrichment analysis
  kegg_results <- gseKEGG(geneList = gene_list_sorted, organism = "hsa", pvalueCutoff = 0.1, pAdjustMethod = "BH")
  if (is.null(kegg_results)) stop("Error: No KEGG enrichment results found")

  # Select top 4 BP, 4 MF, 4 CC terms
  go_df <- data.frame(go_results)
  go_df$ONTOLOGY <- factor(go_df$ONTOLOGY, levels = c("BP", "MF", "CC"))
  go_df_top <- go_df %>% group_by(ONTOLOGY) %>% top_n(-4, p.adjust)

  # KEGG Top 10
  kegg_df <- data.frame(kegg_results) %>% top_n(-10, p.adjust)

  # Create GO Plot
  go_plot <- ggplot(go_df_top, aes(x = reorder(Description, -p.adjust), y = GeneRatio, color = ONTOLOGY)) +
    geom_point(aes(size = -log10(p.adjust))) +
    coord_flip() +
    facet_wrap(~ONTOLOGY, scales = "free") +
    theme_minimal() + labs(title = "GO Enrichment Analysis")

  # Create KEGG Plot
  kegg_plot <- ggplot(kegg_df, aes(x = reorder(Description, -p.adjust), y = GeneRatio)) +
    geom_point(aes(size = -log10(p.adjust), color = Description)) +
    coord_flip() +
    theme_minimal() + labs(title = "KEGG Pathway Enrichment")

  # Combine both plots
  combined_plot <- ggarrange(go_plot, kegg_plot, ncol = 1, nrow = 2)
  ggsave("GO_KEGG_Combined.png", combined_plot, width = 12, height = 10)

  return(list(go_results = go_results, kegg_results = kegg_results, plot = combined_plot))
}



# Example usage:
#plots <- perform_go_kegg_analysis_with_ggplots("LC_SMK.csv")

# Display GO Enrichment Plot (BP, CC, MF)
#print(plots$go_plot)

# Display KEGG Enrichment Plot if available
#if (!is.null(plots$kegg_plot)) {
#  print(plots$kegg_plot)
#}


# Example usage of the function
#results <- perform_go_kegg_analysis_with_plots("covid19.csv")

# To view the GO plot
#print(results$go_plot)

# To view the KEGG plot if available
#if (!is.null(results$kegg_plot)) {
#  print(results$kegg_plot)
#}


