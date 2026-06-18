# Summarize GO Enrichment Results Across Chromosomes
# This script reads GO enrichment results from multiple chromosomes and
# creates a summary showing the top significant hits for each chromosome

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# Set working directory if needed
# setwd("path/to/your/results")

# ============================================================================
# 1. CONFIGURATION
# ============================================================================

# Path to GO results directory (adjust as needed)
# Set to "." if you're already in the directory with the files
go_results_dir <- "."

# Pattern to match GO result files
file_pattern <- ".*_res\\.tab$"

# Significance threshold
bonferroni_threshold <- 0.0001

# Number of top hits to show per chromosome
top_n_hits <- 3

# ============================================================================
# 2. READ ALL GO ENRICHMENT FILES
# ============================================================================

# Find all GO result files
go_files <- list.files(path = go_results_dir, 
                       pattern = file_pattern, 
                       full.names = TRUE)

cat(sprintf("Looking in directory: %s\n", normalizePath(go_results_dir)))
cat(sprintf("Pattern: %s\n", file_pattern))
cat(sprintf("Found %d GO enrichment result files\n\n", length(go_files)))

if (length(go_files) == 0) {
  cat("ERROR: No GO result files found!\n")
  cat("\nTroubleshooting:\n")
  cat("1. Check your current directory:\n")
  cat(sprintf("   Current: %s\n", getwd()))
  cat("2. List files in current directory:\n")
  all_files <- list.files(go_results_dir)
  if (length(all_files) > 0) {
    cat(sprintf("   Found %d total files\n", length(all_files)))
    cat("   First few files:\n")
    print(head(all_files, 10))
  } else {
    cat("   No files found in this directory!\n")
  }
  cat("\n3. Try adjusting the file_pattern variable or go_results_dir path\n")
  stop("No input files found. Please check configuration.")
}

# Read and combine all GO results
all_go_results <- data.frame()

for (go_file in go_files) {
  # Extract chromosome name from filename
  file_basename <- basename(go_file)
  # Pattern: Species_HiC_scaffold_X_res.tab -> extract scaffold X
  chrom <- gsub(".*_([0-9]+)_res\\.tab", "\\1", file_basename)
  
  # Read the file
  tryCatch({
    go_data <- read.table(go_file, 
                          header = TRUE, 
                          sep = "\t",
                          quote = "",
                          comment.char = "",
                          stringsAsFactors = FALSE,
                          fill = TRUE)
    
    # Add chromosome column
    go_data$Chromosome <- paste0("Scaffold_", chrom)
    go_data$Chromosome_Num <- as.numeric(chrom)
    
    # Combine with all results
    all_go_results <- rbind(all_go_results, go_data)
    
    cat(sprintf("  Loaded %s: %d GO terms\n", 
                paste0("Scaffold_", chrom), 
                nrow(go_data)))
    
  }, error = function(e) {
    cat(sprintf("  ERROR reading %s: %s\n", file_basename, e$message))
  })
}

cat(sprintf("\nTotal GO terms loaded: %d\n", nrow(all_go_results)))

if (nrow(all_go_results) == 0) {
  stop("ERROR: No GO data loaded! Check that files are readable and properly formatted.")
}

# Check that required columns exist
required_cols <- c("GO.ID", "Term", "Annotated", "Significant", "bonferroni", "ontology", "FoldChange")
missing_cols <- setdiff(required_cols, colnames(all_go_results))
if (length(missing_cols) > 0) {
  cat("\nERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(colnames(all_go_results), collapse = ", "), "\n")
  stop("Required columns missing from data")
}

# ============================================================================
# 3. FILTER AND SUMMARIZE SIGNIFICANT RESULTS
# ============================================================================

# Filter for significant results (Bonferroni < 0.05)
significant_go <- all_go_results %>%
  filter(bonferroni < bonferroni_threshold) %>%
  arrange(Chromosome_Num, bonferroni)

cat(sprintf("\nSignificant GO terms (Bonferroni < %.2f): %d\n", 
            bonferroni_threshold, 
            nrow(significant_go)))

# Count significant terms per chromosome
sig_per_chrom <- significant_go %>%
  group_by(Chromosome) %>%
  summarise(
    N_Significant = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(N_Significant))

cat("\nSignificant GO terms per chromosome:\n")
print(sig_per_chrom, n = Inf)

# ============================================================================
# 4. GET TOP HITS PER CHROMOSOME
# ============================================================================

top_hits <- significant_go %>%
  group_by(Chromosome) %>%
  slice_min(order_by = bonferroni, n = top_n_hits, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Chromosome_Num, bonferroni) %>%
  select(Chromosome, GO.ID, Term, Annotated, Significant, 
         FoldChange, bonferroni, ontology)

cat(sprintf("\n\n=== TOP %d SIGNIFICANT GO TERMS PER CHROMOSOME ===\n\n", top_n_hits))

# Print formatted output
current_chrom <- ""
for (i in 1:nrow(top_hits)) {
  row <- top_hits[i, ]
  
  # Print chromosome header if it's a new chromosome
  if (row$Chromosome != current_chrom) {
    cat(sprintf("\n%s (%d significant terms total)\n", 
                row$Chromosome,
                sig_per_chrom$N_Significant[sig_per_chrom$Chromosome == row$Chromosome]))
    cat(paste(rep("-", 80), collapse = ""), "\n")
    current_chrom <- row$Chromosome
  }
  
  # Print GO term details
  cat(sprintf("  %d. %s [%s]\n", 
              i - sum(top_hits$Chromosome[1:i-1] != current_chrom) + 1,
              row$Term,
              row$GO.ID))
  cat(sprintf("     Ontology: %s | Genes: %d/%d | Fold Change: %.2fx | Bonferroni P: %.2e\n",
              row$ontology,
              row$Significant,
              row$Annotated,
              row$FoldChange,
              row$bonferroni))
}

# ============================================================================
# 5. CREATE SUMMARY VISUALIZATIONS
# ============================================================================

# --- Figure 1: Number of significant GO terms per chromosome ---

chrom_plot_data <- sig_per_chrom %>%
  mutate(
    Chrom_Num = as.numeric(gsub(".*_([0-9]+).*", "\\1", Chromosome)),
    Chromosome = factor(Chromosome, levels = Chromosome[order(Chrom_Num)])
  )

p1 <- ggplot(chrom_plot_data, aes(x = Chromosome, y = N_Significant)) +
  geom_bar(stat = "identity", fill = "#377EB8") +
  geom_text(aes(label = N_Significant), vjust = -0.5, size = 3) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Significant GO Terms per Chromosome",
    subtitle = sprintf("Bonferroni-corrected P < %.2f", bonferroni_threshold),
    x = "Chromosome",
    y = "Number of Significant GO Terms"
  )

ggsave("go_summary_per_chromosome.pdf", p1, width = 12, height = 6)
ggsave("go_summary_per_chromosome.png", p1, width = 12, height = 6, dpi = 300)
print(p1)

# --- Figure 2: Heatmap variants ---

heatmap_subtitle <- sprintf(
  "Bonferroni < %.2f | >= %d annotated genes",
  bonferroni_threshold,
  min_annotated
)


if (nrow(top_hits) > 0) {
  
  heatmap_data <- top_hits %>%
    mutate(
      Term_Short = ifelse(nchar(Term) > 40, 
                         paste0(substr(Term, 1, 37), "..."), 
                         Term),
      Term_Label = paste0(Term_Short, " (", GO.ID, ")"),
      Chrom_Num = as.numeric(gsub(".*_([0-9]+).*", "\\1", Chromosome)),
      Neg_log10_P = -log10(bonferroni),
      Percent_Significant = 100 * Significant / Annotated
    )
  
  # ==========================================================================
  # VERSION 1: YOUR ORIGINAL IDEA (BEST DEFAULT)
  # Colour = significance, Text = %
  # ==========================================================================
  
  p2_v1 <- ggplot(heatmap_data, aes(x = Chromosome, y = Term_Label, fill = Neg_log10_P)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", Percent_Significant)), size = 2.5) +
    scale_fill_gradient(low = "#FFF7BC", high = "#D95F0E", 
                        name = "-log10(Bonferroni P)") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    ) +
    labs(
      title = sprintf("Top %d GO Terms per Chromosome", top_n_hits),
      subtitle = sprintf(
        "Fill = -log10(Bonferroni P), Text = %% significant genes | Bonferroni < %.2e | >= %d annotated genes",
        bonferroni_threshold,
        min_annotated
      ),
      x = "Chromosome",
      y = "GO Term"
    )
  
  ggsave("go_heatmap_v1_percent_text_pvalue_fill.pdf", p2_v1,
         width = 14, height = max(8, nrow(top_hits) * 0.3))
  print(p2_v1)
  
  
  # ==========================================================================
  # VERSION 2: ADD COUNTS (more informative labels)
  # ==========================================================================
  
  p2_v2 <- ggplot(heatmap_data, aes(x = Chromosome, y = Term_Label, fill = Neg_log10_P)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", 
                                 Percent_Significant, Significant, Annotated)), 
              size = 2.3) +
    scale_fill_gradient(low = "#FFF7BC", high = "#D95F0E",
                        name = "-log10(P)") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    ) +
    labs(
      title = "GO Enrichment Heatmap",
      subtitle = sprintf(
        "Fill = -log10(Bonferroni P), Text = %% significant genes | Bonferroni < %.2e | >= %d annotated genes",
        bonferroni_threshold,
        min_annotated
      ),
      x = "Chromosome",
      y = "GO Term"
    )
  
  ggsave("go_heatmap_v2_percent_counts.pdf", p2_v2,
         width = 14, height = max(8, nrow(top_hits) * 0.3))
  print(p2_v2)
  
  
  # ==========================================================================
  # VERSION 3: SIZE = PERCENT, COLOUR = SIGNIFICANCE
  # (More visual, less text clutter)
  # ==========================================================================
  
  p2_v3 <- ggplot(heatmap_data, aes(x = Chromosome, y = Term_Label)) +
    geom_point(aes(size = Percent_Significant, fill = Neg_log10_P),
               shape = 21, color = "black") +
    scale_fill_gradient(low = "#FFF7BC", high = "#D95F0E",
                        name = "-log10(P)") +
    scale_size(range = c(2, 8), name = "% Significant") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    ) +
    labs(
      title = "GO Enrichment Bubble Plot",
      subtitle = sprintf(
        "Fill = -log10(Bonferroni P), Text = %% significant genes | Bonferroni < %.2e | >= %d annotated genes",
        bonferroni_threshold,
        min_annotated
      ),
      x = "Chromosome",
      y = "GO Term"
    )
  
  ggsave("go_heatmap_v3_bubble.pdf", p2_v3,
         width = 14, height = max(8, nrow(top_hits) * 0.3))
  print(p2_v3)
  
  
  # ==========================================================================
  # VERSION 4: FACET BY ONTOLOGY (BP / MF / CC)
  # ==========================================================================
  
  p2_v4 <- ggplot(heatmap_data, aes(x = Chromosome, y = Term_Label, fill = Neg_log10_P)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", Percent_Significant)), size = 2.3) +
    scale_fill_gradient(low = "#FFF7BC", high = "#D95F0E",
                        name = "-log10(P)") +
    facet_wrap(~ontology, scales = "free_y") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(
      title = "GO Enrichment by Ontology",
      subtitle = sprintf(
        "Fill = -log10(Bonferroni P), Text = %% significant genes | Bonferroni < %.2e | >= %d annotated genes",
        bonferroni_threshold,
        min_annotated
      ),
      x = "Chromosome",
      y = "GO Term"
    )
  
  ggsave("go_heatmap_v4_faceted.pdf", p2_v4,
         width = 14, height = max(8, nrow(top_hits) * 0.35))
  print(p2_v4)
  
  
  # ==========================================================================
  # VERSION 5: FILTER LOW ANNOTATION (more robust)
  # ==========================================================================
  
  heatmap_filtered <- heatmap_data %>%
    filter(Annotated >= 10)
  
  p2_v5 <- ggplot(heatmap_filtered, aes(x = Chromosome, y = Term_Label, fill = Neg_log10_P)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", Percent_Significant)), size = 2.5) +
    scale_fill_gradient(low = "#FFF7BC", high = "#D95F0E",
                        name = "-log10(P)") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    ) +
    labs(
      title = "GO Enrichment (Filtered)",
      subtitle = sprintf(
        "Fill = -log10(Bonferroni P), Text = %% significant genes | Bonferroni < %.2e | >= %d annotated genes",
        bonferroni_threshold,
        min_annotated
      ),
      x = "Chromosome",
      y = "GO Term"
    )
  
  ggsave("go_heatmap_v5_filtered.pdf", p2_v5,
         width = 14, height = max(8, nrow(heatmap_filtered) * 0.3))
  print(p2_v5)
}


# --- Figure 3: Distribution of ontologies ---

ontology_summary <- significant_go %>%
  group_by(Chromosome, ontology) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(
    Chrom_Num = as.numeric(gsub(".*_([0-9]+).*", "\\1", Chromosome)),
    Chromosome = factor(Chromosome, levels = unique(Chromosome[order(Chrom_Num)]))
  )

p3 <- ggplot(ontology_summary, aes(x = Chromosome, y = Count, fill = ontology)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c("BP" = "#E41A1C", "MF" = "#377EB8", "CC" = "#4DAF4A"),
    labels = c("BP" = "Biological Process", 
               "MF" = "Molecular Function", 
               "CC" = "Cellular Component")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title = "GO Term Ontology Distribution per Chromosome",
    subtitle = sprintf("Bonferroni-corrected P < %.2f", bonferroni_threshold),
    x = "Chromosome",
    y = "Number of Significant GO Terms",
    fill = "Ontology"
  )

ggsave("go_summary_ontology.pdf", p3, width = 12, height = 6)
ggsave("go_summary_ontology.png", p3, width = 12, height = 6, dpi = 300)
print(p3)

# ============================================================================
# 6. EXPORT RESULTS
# ============================================================================

# Export full significant results
write.csv(significant_go, "go_significant_all_chromosomes.csv", row.names = FALSE)

# Export top hits summary
write.csv(top_hits, "go_top_hits_per_chromosome.csv", row.names = FALSE)

# Export summary statistics
summary_stats <- all_go_results %>%
  group_by(Chromosome) %>%
  summarise(
    Total_GO_Terms = n(),
    Significant_Terms = sum(bonferroni < bonferroni_threshold),
    Pct_Significant = 100 * sum(bonferroni < bonferroni_threshold) / n(),
    Min_Bonferroni = min(bonferroni),
    Median_FoldChange = median(FoldChange[bonferroni < bonferroni_threshold], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Significant_Terms))

write.csv(summary_stats, "go_summary_statistics.csv", row.names = FALSE)

# Create a detailed report
report_file <- "go_enrichment_summary_report.txt"
sink(report_file)

cat("=" , rep("=", 78), "\n", sep = "")
cat("GO ENRICHMENT SUMMARY REPORT\n")
cat("=" , rep("=", 78), "\n\n", sep = "")

cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results Directory:", go_results_dir, "\n")
cat("Bonferroni Threshold:", bonferroni_threshold, "\n")
cat("Top Hits per Chromosome:", top_n_hits, "\n\n")

cat("OVERALL SUMMARY\n")
cat("-" , rep("-", 78), "\n", sep = "")
cat("Total chromosomes analyzed:", length(unique(all_go_results$Chromosome)), "\n")
cat("Total GO terms tested:", nrow(all_go_results), "\n")
cat("Significant GO terms:", nrow(significant_go), "\n")
cat("Percentage significant:", 
    sprintf("%.1f%%", 100 * nrow(significant_go) / nrow(all_go_results)), "\n\n")

cat("CHROMOSOMES WITH MOST ENRICHMENT\n")
cat("-" , rep("-", 78), "\n", sep = "")
print(summary_stats, n = Inf)

cat("\n\nTOP GO TERMS PER CHROMOSOME\n")
cat("-" , rep("-", 78), "\n", sep = "")

current_chrom <- ""
for (i in 1:nrow(top_hits)) {
  row <- top_hits[i, ]
  
  if (row$Chromosome != current_chrom) {
    cat("\n", row$Chromosome, "\n", sep = "")
    cat(paste(rep("-", 78), collapse = ""), "\n")
    current_chrom <- row$Chromosome
  }
  
  cat(sprintf("%d. %s [%s]\n", 
              i - sum(top_hits$Chromosome[1:i-1] != current_chrom) + 1,
              row$Term,
              row$GO.ID))
  cat(sprintf("   Ontology: %s | Enrichment: %.2fx | P: %.2e\n",
              row$ontology,
              row$FoldChange,
              row$bonferroni))
}

# ============================================================================
# 7. CREATE FINAL SUMMARY TABLES
# ============================================================================

cat("\nGenerating final summary tables...\n")

# ---- 7.1 Build master table ----

final_summary <- significant_go %>%
  mutate(
    Percent_Significant = 100 * Significant / Annotated,
    Term_Clean = paste0(Term, " (", GO.ID, ")")
  ) %>%
  # Join number of significant GO terms per chromosome
  left_join(sig_per_chrom, by = "Chromosome") %>%
  select(
    Chromosome,
    Chromosome_Num,
    Term_Clean,
    GO.ID,
    Term,
    ontology,
    Significant,
    Annotated,
    Percent_Significant,
    FoldChange,
    bonferroni,
    N_Significant
  )

# ---- 7.2 Robust ordering ----
# Numeric chromosomes first, then non-numeric

final_summary <- final_summary %>%
  arrange(
    is.na(Chromosome_Num),
    Chromosome_Num,
    bonferroni
  )

# ---- 7.3 Write full table ----

write.csv(final_summary, 
          "go_final_summary_table_full.csv", 
          row.names = FALSE)

cat("  - Full summary table written: go_final_summary_table_full.csv\n")

# ---- 7.4 Sorted versions ----

# By highest % significant (global)
final_summary_pct <- final_summary %>%
  arrange(desc(Percent_Significant), bonferroni)

write.csv(final_summary_pct, 
          "go_final_summary_sorted_by_percent.csv", 
          row.names = FALSE)

# By chromosome, then % significant
final_summary_by_chr <- final_summary %>%
  group_by(Chromosome) %>%
  arrange(desc(Percent_Significant), .by_group = TRUE) %>%
  ungroup()

write.csv(final_summary_by_chr, 
          "go_final_summary_by_chromosome.csv", 
          row.names = FALSE)

cat("  - Sorted tables written\n")

# ---- 7.5 Human-readable text output ----

final_summary_readable <- final_summary %>%
  mutate(
    Summary_Line = sprintf(
      "%s\t%.1f%% (%d/%d)\t%s\tP=%.2e\tFC=%.2fx",
      Term_Clean,
      Percent_Significant,
      Significant,
      Annotated,
      Chromosome,
      bonferroni,
      FoldChange
    )
  )

writeLines(final_summary_readable$Summary_Line,
           "go_final_summary_readable.txt")

cat("  - Readable summary written: go_final_summary_readable.txt\n")

# ---- 7.6 Optional: filter low annotation terms ----

min_annotated <- 10

final_summary_filtered <- final_summary %>%
  filter(Annotated >= min_annotated)

write.csv(final_summary_filtered,
          "go_final_summary_filtered_annotated_ge10.csv",
          row.names = FALSE)

cat(sprintf("  - Filtered table (Annotated >= %d) written\n", min_annotated))

# ---- 7.7 Optional: top terms per chromosome (by %) ----

top_terms_per_chr <- final_summary %>%
  group_by(Chromosome) %>%
  slice_max(order_by = Percent_Significant, n = top_n_hits, with_ties = FALSE) %>%
  ungroup()

write.csv(top_terms_per_chr,
          "go_top_terms_by_percent_per_chromosome.csv",
          row.names = FALSE)

cat("  - Top terms per chromosome (by %) written\n")

# ---- 7.8 Summary message ----

cat("\nFinal summary tables complete.\n\n")

sink()

cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("\nGenerated files:\n")
cat("  FIGURES:\n")
cat("    - go_summary_per_chromosome.pdf/png (bar chart of significant terms)\n")
cat("    - go_summary_heatmap.pdf/png (heatmap of top terms)\n")
cat("    - go_summary_ontology.pdf/png (ontology distribution)\n")
cat("  DATA:\n")
cat("    - go_significant_all_chromosomes.csv (all significant GO terms)\n")
cat("    - go_top_hits_per_chromosome.csv (top hits summary)\n")
cat("    - go_summary_statistics.csv (summary statistics per chromosome)\n")
cat("    - go_enrichment_summary_report.txt (formatted text report)\n")
cat("\n")
