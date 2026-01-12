# Orthogroup Analysis: Queen vs Worker Up-regulated Genes Across Species
# This script analyzes differential gene expression across three species
# using orthogroup information to identify conserved regulatory patterns

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(VennDiagram)
  library(grid)
  library(futile.logger)
})

# Set working directory (adjust as needed)
# setwd("your/path/here")

# ============================================================================
# 1. READ DATA
# ============================================================================

# Read orthogroups file
orthogroups <- read.delim("Orthogroups.tsv", 
                          header = TRUE, 
                          stringsAsFactors = FALSE)

# Define species abbreviations (adjust these to match your filenames)
species <- c("BJ", "LF", "PD")  # Belonogaster juncea, Liostenogaster flavolineata, Polistes dominula
species_names <- c("B. juncea", "L. flavolineata", "P. dominula")

# Read all gene lists
gene_lists <- list()

for (sp in species) {
  # Queen-up genes
  queen_file <- paste0(sp, "_QUEEN_UP.txt")
  if (file.exists(queen_file)) {
    gene_lists[[paste0(sp, "_QUEEN_UP")]] <- read.table(queen_file, 
                                                         header = FALSE, 
                                                         stringsAsFactors = FALSE)$V1
  }
  
  # Worker-up genes
  worker_file <- paste0(sp, "_WORKER_UP.txt")
  if (file.exists(worker_file)) {
    gene_lists[[paste0(sp, "_WORKER_UP")]] <- read.table(worker_file, 
                                                          header = FALSE, 
                                                          stringsAsFactors = FALSE)$V1
  }
}

# Print summary
cat("Gene lists loaded:\n")
for (name in names(gene_lists)) {
  cat(sprintf("  %s: %d genes\n", name, length(gene_lists[[name]])))
}

# ============================================================================
# 2. PROCESS ORTHOGROUPS
# ============================================================================

# Function to extract gene IDs from orthogroup entries (remove isoform suffixes)
extract_gene_ids <- function(gene_string) {
  if (is.na(gene_string) || gene_string == "") return(character(0))
  genes <- strsplit(gene_string, ", ")[[1]]
  # Remove isoform suffixes (.t1, .t2, etc.)
  genes <- gsub("\\.t[0-9]+$", "", genes)
  unique(genes)
}

# Process orthogroups to create a mapping
orthogroup_data <- orthogroups %>%
  rowwise() %>%
  mutate(
    BJ_genes = list(extract_gene_ids(Belonogaster_juncea)),
    LF_genes = list(extract_gene_ids(Liostenogaster_flavolineata)),
    PD_genes = list(extract_gene_ids(Polistes_dominula))
  ) %>%
  ungroup()

# ============================================================================
# 3. MAP GENES TO ORTHOGROUPS
# ============================================================================

# Function to find orthogroups for a given gene list and species
map_genes_to_orthogroups <- function(gene_list, species_col, ortho_data) {
  orthogroups_found <- c()
  
  for (og_idx in 1:nrow(ortho_data)) {
    genes_in_og <- ortho_data[[species_col]][[og_idx]]
    if (any(gene_list %in% genes_in_og)) {
      orthogroups_found <- c(orthogroups_found, ortho_data$Orthogroup[og_idx])
    }
  }
  
  unique(orthogroups_found)
}

# Map all gene lists to orthogroups
orthogroup_assignments <- list()

for (sp_idx in 1:length(species)) {
  sp <- species[sp_idx]
  species_col <- paste0(sp, "_genes")
  
  # Queen-up
  queen_key <- paste0(sp, "_QUEEN_UP")
  if (queen_key %in% names(gene_lists)) {
    orthogroup_assignments[[queen_key]] <- 
      map_genes_to_orthogroups(gene_lists[[queen_key]], species_col, orthogroup_data)
  }
  
  # Worker-up
  worker_key <- paste0(sp, "_WORKER_UP")
  if (worker_key %in% names(gene_lists)) {
    orthogroup_assignments[[worker_key]] <- 
      map_genes_to_orthogroups(gene_lists[[worker_key]], species_col, orthogroup_data)
  }
}

# ============================================================================
# 4. CREATE PRESENCE/ABSENCE MATRIX
# ============================================================================

# Get all unique orthogroups that appear in at least one condition
all_orthogroups <- unique(unlist(orthogroup_assignments))

# Create matrix: rows = orthogroups, columns = species x condition
matrix_data <- matrix(0, 
                      nrow = length(all_orthogroups), 
                      ncol = length(orthogroup_assignments))

rownames(matrix_data) <- all_orthogroups
colnames(matrix_data) <- names(orthogroup_assignments)

# Fill in the matrix
for (col_idx in 1:length(orthogroup_assignments)) {
  col_name <- names(orthogroup_assignments)[col_idx]
  ogs <- orthogroup_assignments[[col_name]]
  matrix_data[ogs, col_idx] <- 1
}

# ============================================================================
# 5. ANALYZE PATTERNS
# ============================================================================

# Count orthogroups by pattern
pattern_summary <- data.frame(
  Orthogroup = all_orthogroups,
  matrix_data,
  check.names = FALSE
) %>%
  mutate(
    Total_species = rowSums(select(., contains("_QUEEN_UP")) | select(., contains("_WORKER_UP"))),
    Queen_count = rowSums(select(., contains("_QUEEN_UP"))),
    Worker_count = rowSums(select(., contains("_WORKER_UP"))),
    Pattern = case_when(
      Queen_count == 3 ~ "Queen-up (all 3 species)",
      Worker_count == 3 ~ "Worker-up (all 3 species)",
      Queen_count == 2 ~ "Queen-up (2 species)",
      Worker_count == 2 ~ "Worker-up (2 species)",
      Queen_count > 0 & Worker_count > 0 ~ "Mixed (Queen & Worker)",
      TRUE ~ "Single species"
    )
  )

# Print summary statistics
cat("\n\n=== PATTERN SUMMARY ===\n")
pattern_counts <- table(pattern_summary$Pattern)
print(pattern_counts)

# Identify conserved orthogroups (present in all 3 species)
conserved_queen <- pattern_summary %>%
  filter(Queen_count == 3, Worker_count == 0) %>%
  pull(Orthogroup)

conserved_worker <- pattern_summary %>%
  filter(Worker_count == 3, Queen_count == 0) %>%
  pull(Orthogroup)

mixed_regulation <- pattern_summary %>%
  filter(Queen_count > 0, Worker_count > 0) %>%
  pull(Orthogroup)

cat(sprintf("\nConserved Queen-up orthogroups: %d\n", length(conserved_queen)))
cat(sprintf("Conserved Worker-up orthogroups: %d\n", length(conserved_worker)))
cat(sprintf("Mixed regulation orthogroups: %d\n", length(mixed_regulation)))

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================

# ============================================================================
# 6a. STATISTICAL TESTING OF OVERLAPS
# ============================================================================

# Function to perform hypergeometric test for overlap significance
test_overlap_significance <- function(set1, set2, universe_size) {
  # set1, set2: vectors of orthogroup IDs
  # universe_size: total number of orthogroups
  
  overlap <- length(intersect(set1, set2))
  size1 <- length(set1)
  size2 <- length(set2)
  
  # Hypergeometric test
  # P(X >= overlap) where X ~ Hypergeometric(universe_size, size1, size2)
  p_value <- phyper(overlap - 1, size1, universe_size - size1, size2, lower.tail = FALSE)
  
  # Calculate expected overlap under null hypothesis
  expected <- (size1 * size2) / universe_size
  
  # Fold enrichment
  fold_enrichment <- ifelse(expected > 0, overlap / expected, NA)
  
  return(list(
    overlap = overlap,
    expected = expected,
    fold_enrichment = fold_enrichment,
    p_value = p_value,
    size1 = size1,
    size2 = size2
  ))
}

# Get universe size (total unique orthogroups)
universe_size <- nrow(orthogroup_data)

# Test all pairwise comparisons within and between species
overlap_tests <- list()

# Within-species tests (Queen vs Worker)
cat("\n\n=== WITHIN-SPECIES OVERLAP TESTS ===\n")
for (i in 1:length(species)) {
  sp <- species[i]
  queen_key <- paste0(sp, "_QUEEN_UP")
  worker_key <- paste0(sp, "_WORKER_UP")
  
  if (queen_key %in% names(orthogroup_assignments) && 
      worker_key %in% names(orthogroup_assignments)) {
    
    test_result <- test_overlap_significance(
      orthogroup_assignments[[queen_key]],
      orthogroup_assignments[[worker_key]],
      universe_size
    )
    
    overlap_tests[[paste(sp, "Queen_vs_Worker")]] <- test_result
    
    cat(sprintf("\n%s (Queen vs Worker):\n", species_names[i]))
    cat(sprintf("  Observed overlap: %d\n", test_result$overlap))
    cat(sprintf("  Expected overlap: %.1f\n", test_result$expected))
    cat(sprintf("  Fold enrichment: %.2f\n", test_result$fold_enrichment))
    cat(sprintf("  P-value: %.2e %s\n", test_result$p_value, 
                ifelse(test_result$p_value < 0.001, "***",
                       ifelse(test_result$p_value < 0.01, "**",
                              ifelse(test_result$p_value < 0.05, "*", "ns")))))
  }
}

# Between-species tests (same caste)
cat("\n\n=== BETWEEN-SPECIES OVERLAP TESTS ===\n")
for (caste in c("QUEEN", "WORKER")) {
  cat(sprintf("\n%s-biased genes:\n", caste))
  
  for (i in 1:(length(species)-1)) {
    for (j in (i+1):length(species)) {
      sp1 <- species[i]
      sp2 <- species[j]
      key1 <- paste0(sp1, "_", caste, "_UP")
      key2 <- paste0(sp2, "_", caste, "_UP")
      
      if (key1 %in% names(orthogroup_assignments) && 
          key2 %in% names(orthogroup_assignments)) {
        
        test_result <- test_overlap_significance(
          orthogroup_assignments[[key1]],
          orthogroup_assignments[[key2]],
          universe_size
        )
        
        test_name <- paste(sp1, "vs", sp2, caste)
        overlap_tests[[test_name]] <- test_result
        
        cat(sprintf("\n  %s vs %s:\n", species_names[i], species_names[j]))
        cat(sprintf("    Observed overlap: %d\n", test_result$overlap))
        cat(sprintf("    Expected overlap: %.1f\n", test_result$expected))
        cat(sprintf("    Fold enrichment: %.2f\n", test_result$fold_enrichment))
        cat(sprintf("    P-value: %.2e %s\n", test_result$p_value, 
                    ifelse(test_result$p_value < 0.001, "***",
                           ifelse(test_result$p_value < 0.01, "**",
                                  ifelse(test_result$p_value < 0.05, "*", "ns")))))
      }
    }
  }
}

# Three-way overlaps
cat("\n\n=== THREE-WAY OVERLAP TESTS ===\n")
for (caste in c("QUEEN", "WORKER")) {
  keys <- paste0(species, "_", caste, "_UP")
  
  if (all(keys %in% names(orthogroup_assignments))) {
    # Get pairwise intersections first
    int_12 <- intersect(orthogroup_assignments[[keys[1]]], 
                        orthogroup_assignments[[keys[2]]])
    three_way <- intersect(int_12, orthogroup_assignments[[keys[3]]])
    
    # Test if three-way overlap is significant given the pairwise overlaps
    # Use the smallest pairwise intersection as the "universe"
    pairwise_sizes <- c(
      length(intersect(orthogroup_assignments[[keys[1]]], orthogroup_assignments[[keys[2]]])),
      length(intersect(orthogroup_assignments[[keys[1]]], orthogroup_assignments[[keys[3]]])),
      length(intersect(orthogroup_assignments[[keys[2]]], orthogroup_assignments[[keys[3]]]))
    )
    
    cat(sprintf("\n%s-biased (all 3 species):\n", caste))
    cat(sprintf("  Three-way overlap: %d orthogroups\n", length(three_way)))
    cat(sprintf("  Pairwise overlaps: %s\n", paste(pairwise_sizes, collapse = ", ")))
  }
}

# Convert overlap tests to data frame for export
overlap_test_df <- do.call(rbind, lapply(names(overlap_tests), function(name) {
  test <- overlap_tests[[name]]
  data.frame(
    Comparison = name,
    Observed = test$overlap,
    Expected = round(test$expected, 2),
    Fold_Enrichment = round(test$fold_enrichment, 2),
    P_value = test$p_value,
    Significance = ifelse(test$p_value < 0.001, "***",
                         ifelse(test$p_value < 0.01, "**",
                                ifelse(test$p_value < 0.05, "*", "ns"))),
    stringsAsFactors = FALSE
  )
}))

write.csv(overlap_test_df, "overlap_significance_tests.csv", row.names = FALSE)

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================

# --- Figure 1: UpSet Plot (showing intersections) ---
# Using a manual approach since UpSetR might not be available

# Calculate intersections for each caste separately
calculate_intersections <- function(caste) {
  sets <- list()
  for (sp in species) {
    key <- paste0(sp, "_", caste, "_UP")
    if (key %in% names(orthogroup_assignments)) {
      sets[[sp]] <- orthogroup_assignments[[key]]
    }
  }
  
  # Calculate all possible intersections
  intersections <- list()
  
  # Single species
  for (sp in species) {
    if (sp %in% names(sets)) {
      exclusive <- sets[[sp]]
      for (other_sp in setdiff(species, sp)) {
        if (other_sp %in% names(sets)) {
          exclusive <- setdiff(exclusive, sets[[other_sp]])
        }
      }
      intersections[[sp]] <- exclusive
    }
  }
  
  # Pairwise intersections
  for (i in 1:(length(species)-1)) {
    for (j in (i+1):length(species)) {
      sp1 <- species[i]
      sp2 <- species[j]
      if (sp1 %in% names(sets) && sp2 %in% names(sets)) {
        pair_int <- intersect(sets[[sp1]], sets[[sp2]])
        # Exclude those also in third species
        other_sp <- setdiff(species, c(sp1, sp2))
        for (osp in other_sp) {
          if (osp %in% names(sets)) {
            pair_int <- setdiff(pair_int, sets[[osp]])
          }
        }
        intersections[[paste(sp1, sp2, sep = " & ")]] <- pair_int
      }
    }
  }
  
  # Three-way intersection
  if (all(species %in% names(sets))) {
    three_way <- Reduce(intersect, sets)
    intersections[[paste(species, collapse = " & ")]] <- three_way
  }
  
  return(intersections)
}

queen_intersections <- calculate_intersections("QUEEN")
worker_intersections <- calculate_intersections("WORKER")

# Create intersection plot data
create_intersection_df <- function(intersections, caste) {
  df <- data.frame(
    Intersection = names(intersections),
    Count = sapply(intersections, length),
    Caste = caste,
    stringsAsFactors = FALSE
  )
  df$Intersection <- factor(df$Intersection, levels = df$Intersection[order(df$Count, decreasing = TRUE)])
  return(df)
}

queen_int_df <- create_intersection_df(queen_intersections, "Queen")
worker_int_df <- create_intersection_df(worker_intersections, "Worker")
all_int_df <- rbind(queen_int_df, worker_int_df)

# Add significance annotations to the plot data
# Match intersection names to test results
all_int_df$Significance <- ""
for (i in 1:nrow(all_int_df)) {
  int_name <- as.character(all_int_df$Intersection[i])
  caste <- all_int_df$Caste[i]
  
  # For two-way intersections
  if (grepl(" & ", int_name) && !grepl("&.*&", int_name)) {
    species_pair <- strsplit(int_name, " & ")[[1]]
    
    # Find matching test
    test_name <- paste(species_pair[1], "vs", species_pair[2], toupper(caste))
    if (test_name %in% names(overlap_tests)) {
      p_val <- overlap_tests[[test_name]]$p_value
      all_int_df$Significance[i] <- ifelse(p_val < 0.001, "***",
                                           ifelse(p_val < 0.01, "**",
                                                  ifelse(p_val < 0.05, "*", "")))
    }
  }
}

# Plot intersections with significance
p1 <- ggplot(all_int_df, aes(x = Intersection, y = Count, fill = Caste)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Significance), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Queen" = "#E41A1C", "Worker" = "#377EB8")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Orthogroup Intersections Across Species",
    subtitle = "Number of orthogroups shared between species for each caste\n(*, **, *** indicate P < 0.05, 0.01, 0.001 respectively)",
    x = "Species Combination",
    y = "Number of Orthogroups",
    fill = "Caste Bias"
  )

ggsave("figure1_intersections.pdf", p1, width = 10, height = 6)
ggsave("figure1_intersections.png", p1, width = 10, height = 6, dpi = 300)
print(p1)

# --- Figure 2: Heatmap of Orthogroup Presence ---

# Filter to orthogroups present in at least 2 conditions
filtered_matrix <- matrix_data[rowSums(matrix_data) >= 2, ]

# Reorder columns for better visualization
col_order <- c(
  paste0(species, "_QUEEN_UP"),
  paste0(species, "_WORKER_UP")
)
col_order <- col_order[col_order %in% colnames(filtered_matrix)]
filtered_matrix <- filtered_matrix[, col_order]

# Create nicer column names
nice_names <- gsub("_", " ", colnames(filtered_matrix))
nice_names <- gsub("BJ", "B. juncea", nice_names)
nice_names <- gsub("LF", "L. flavolineata", nice_names)
nice_names <- gsub("PD", "P. dominula", nice_names)

# Create annotations
col_annotation <- data.frame(
  Species = factor(rep(species_names, 2), levels = species_names),
  Caste = factor(rep(c("Queen", "Worker"), each = length(species)), levels = c("Queen", "Worker"))
)
rownames(col_annotation) <- colnames(filtered_matrix)

# Create color scheme
col_colors <- list(
  Species = setNames(brewer.pal(length(species_names), "Set2"), species_names),
  Caste = c("Queen" = "#E41A1C", "Worker" = "#377EB8")
)

# Create heatmap
pdf("figure2_heatmap.pdf", width = 8, height = 10)
Heatmap(
  filtered_matrix,
  name = "Present",
  col = c("0" = "white", "1" = "black"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_labels = nice_names,
  column_names_rot = 45,
  column_names_side = "bottom",
  top_annotation = HeatmapAnnotation(
    df = col_annotation,
    col = col_colors,
    show_legend = TRUE
  ),
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("Absent", "Present")
  ),
  column_title = "Orthogroup Presence Across Species and Castes",
  row_title = sprintf("Orthogroups (n=%d)", nrow(filtered_matrix))
)
dev.off()

# --- Figure 3: Venn Diagram Style Representation ---

# Create a summary plot showing conservation patterns
pattern_plot_data <- pattern_summary %>%
  count(Pattern) %>%
  mutate(
    Pattern = factor(Pattern, levels = Pattern[order(n, decreasing = TRUE)]),
    Color = case_when(
      grepl("Queen.*all 3", Pattern) ~ "#E41A1C",
      grepl("Worker.*all 3", Pattern) ~ "#377EB8",
      grepl("Mixed", Pattern) ~ "#984EA3",
      grepl("Queen.*2", Pattern) ~ "#FF7F00",
      grepl("Worker.*2", Pattern) ~ "#4DAF4A",
      TRUE ~ "#999999"
    )
  )

p3 <- ggplot(pattern_plot_data, aes(x = Pattern, y = n, fill = Pattern)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = setNames(pattern_plot_data$Color, pattern_plot_data$Pattern)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Conservation Patterns of Differentially Expressed Orthogroups",
    x = "Pattern",
    y = "Number of Orthogroups"
  )

ggsave("figure3_patterns.pdf", p3, width = 10, height = 6)
ggsave("figure3_patterns.png", p3, width = 10, height = 6, dpi = 300)
print(p3)

# --- Figure 4: Comparison Matrix ---

# Create a matrix showing overlap between queen and worker for each species
comparison_data <- data.frame(
  Species = rep(species_names, 3),
  Category = rep(c("Queen only", "Both", "Worker only"), each = length(species)),
  Count = 0
)

for (i in 1:length(species)) {
  sp <- species[i]
  queen_key <- paste0(sp, "_QUEEN_UP")
  worker_key <- paste0(sp, "_WORKER_UP")
  
  if (queen_key %in% names(orthogroup_assignments) && 
      worker_key %in% names(orthogroup_assignments)) {
    
    queen_ogs <- orthogroup_assignments[[queen_key]]
    worker_ogs <- orthogroup_assignments[[worker_key]]
    
    queen_only <- setdiff(queen_ogs, worker_ogs)
    worker_only <- setdiff(worker_ogs, queen_ogs)
    both <- intersect(queen_ogs, worker_ogs)
    
    comparison_data$Count[comparison_data$Species == species_names[i] & 
                           comparison_data$Category == "Queen only"] <- length(queen_only)
    comparison_data$Count[comparison_data$Species == species_names[i] & 
                           comparison_data$Category == "Both"] <- length(both)
    comparison_data$Count[comparison_data$Species == species_names[i] & 
                           comparison_data$Category == "Worker only"] <- length(worker_only)
  }
}

comparison_data$Category <- factor(comparison_data$Category, 
                                   levels = c("Queen only", "Both", "Worker only"))

p4 <- ggplot(comparison_data, aes(x = Species, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Queen only" = "#E41A1C", 
                                "Both" = "#984EA3", 
                                "Worker only" = "#377EB8")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "Queen vs Worker Orthogroup Overlap Within Each Species",
    x = "Species",
    y = "Number of Orthogroups",
    fill = "Regulation Pattern"
  )

ggsave("figure4_within_species.pdf", p4, width = 8, height = 6)
ggsave("figure4_within_species.png", p4, width = 8, height = 6, dpi = 300)
print(p4)

# --- Figure 5: Overlap Significance Plot ---

# Prepare data for visualization
sig_plot_data <- overlap_test_df %>%
  mutate(
    Comparison_Type = case_when(
      grepl("Queen_vs_Worker", Comparison) ~ "Within-species\n(Queen vs Worker)",
      grepl("QUEEN", Comparison) ~ "Between-species\n(Queen-biased)",
      grepl("WORKER", Comparison) ~ "Between-species\n(Worker-biased)",
      TRUE ~ "Other"
    ),
    Comparison_Label = gsub("_", " ", Comparison),
    Comparison_Label = gsub(" QUEEN| WORKER", "", Comparison_Label),
    Significant = P_value < 0.05,
    Neg_log10_P = -log10(P_value)
  ) %>%
  filter(Comparison_Type != "Other")

# Create fold enrichment plot
p5 <- ggplot(sig_plot_data, aes(x = reorder(Comparison_Label, Fold_Enrichment), 
                                 y = Fold_Enrichment, 
                                 fill = Significant)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = Significance), hjust = -0.2, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#999999"),
                    labels = c("TRUE" = "P < 0.05", "FALSE" = "P ≥ 0.05")) +
  facet_wrap(~ Comparison_Type, scales = "free_y", ncol = 1) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "Overlap Enrichment Analysis",
    subtitle = "Fold enrichment of observed vs expected orthogroup overlaps\n(dashed line = expected overlap by chance)",
    x = "Comparison",
    y = "Fold Enrichment\n(Observed / Expected)",
    fill = "Significance"
  )

ggsave("figure5_enrichment.pdf", p5, width = 10, height = 10)
ggsave("figure5_enrichment.png", p5, width = 10, height = 10, dpi = 300)
print(p5)

# Create a volcano-style plot (fold enrichment vs significance)
p6 <- ggplot(sig_plot_data, aes(x = Fold_Enrichment, y = Neg_log10_P, 
                                 color = Comparison_Type, 
                                 size = Observed)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_text(aes(label = ifelse(Significant, Comparison_Label, "")), 
            hjust = 0, vjust = 0, size = 3, show.legend = FALSE, nudge_x = 0.1) +
  scale_color_brewer(palette = "Set1") +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right") +
  labs(
    title = "Overlap Enrichment: Magnitude vs Significance",
    subtitle = "Dashed lines indicate fold enrichment = 1 and P = 0.05",
    x = "Fold Enrichment (Observed / Expected)",
    y = "-log10(P-value)",
    color = "Comparison Type",
    size = "Observed\nOverlap"
  )

ggsave("figure6_volcano.pdf", p6, width = 12, height = 8)
ggsave("figure6_volcano.png", p6, width = 12, height = 8, dpi = 300)
print(p6)

# --- Figure 7: Venn Diagrams with Statistical Annotations ---

# Load VennDiagram package
suppressPackageStartupMessages(library(VennDiagram))
library(futile.logger)

# Suppress VennDiagram log messages
flog.threshold(ERROR)

# Prepare data for Venn diagrams
queen_sets <- list()
worker_sets <- list()

for (i in 1:length(species)) {
  sp <- species[i]
  sp_name <- species_names[i]
  
  queen_key <- paste0(sp, "_QUEEN_UP")
  worker_key <- paste0(sp, "_WORKER_UP")
  
  if (queen_key %in% names(orthogroup_assignments)) {
    queen_sets[[sp_name]] <- orthogroup_assignments[[queen_key]]
  }
  
  if (worker_key %in% names(orthogroup_assignments)) {
    worker_sets[[sp_name]] <- orthogroup_assignments[[worker_key]]
  }
}

# Function to create annotated Venn diagram with stats
create_venn_with_stats <- function(sets, title, colors, test_results) {
  # Create the base Venn diagram
  venn_plot <- venn.diagram(
    x = sets,
    category.names = names(sets),
    filename = NULL,
    output = TRUE,
    
    # Output features
    imagetype = "png",
    height = 2000,
    width = 2000,
    resolution = 300,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = colors,
    alpha = 0.5,
    
    # Numbers
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
    cat.dist = c(0.055, 0.055, 0.055),
    
    # Title
    main = title,
    main.cex = 1.8,
    main.fontface = "bold",
    main.fontfamily = "sans",
    main.pos = c(0.5, 1.05)
  )
  
  return(venn_plot)
}

# Create Queen-upregulated Venn diagram
venn_queen <- create_venn_with_stats(
  queen_sets,
  "Queen-upregulated Orthogroups",
  c("#E41A1C", "#FF7F00", "#FFFF33"),
  overlap_tests
)

# Save Queen Venn diagram
pdf("figure7_venn_queen.pdf", width = 10, height = 8)
grid.newpage()
grid.draw(venn_queen)
dev.off()

png("figure7_venn_queen.png", width = 2400, height = 1920, res = 300)
grid.newpage()
grid.draw(venn_queen)
dev.off()

# Create Worker-upregulated Venn diagram
venn_worker <- create_venn_with_stats(
  worker_sets,
  "Worker-upregulated Orthogroups",
  c("#377EB8", "#4DAF4A", "#984EA3"),
  overlap_tests
)

# Save Worker Venn diagram
pdf("figure7_venn_worker.pdf", width = 10, height = 8)
grid.newpage()
grid.draw(venn_worker)
dev.off()

png("figure7_venn_worker.png", width = 2400, height = 1920, res = 300)
grid.newpage()
grid.draw(venn_worker)
dev.off()

# Print summary of overlaps with statistics
cat("\n\n=== VENN DIAGRAM OVERLAPS WITH STATISTICS ===\n")
cat("\nQueen-upregulated orthogroups:\n")
for (name in names(queen_sets)) {
  cat(sprintf("  %s: %d orthogroups\n", name, length(queen_sets[[name]])))
}

# Calculate pairwise overlaps for queens
if (length(queen_sets) >= 2) {
  sp_names <- names(queen_sets)
  for (i in 1:(length(sp_names)-1)) {
    for (j in (i+1):length(sp_names)) {
      overlap <- length(intersect(queen_sets[[sp_names[i]]], queen_sets[[sp_names[j]]]))
      test_name <- paste(species[i], "vs", species[j], "QUEEN")
      if (test_name %in% names(overlap_tests)) {
        test <- overlap_tests[[test_name]]
        cat(sprintf("  %s & %s: %d orthogroups (%.2fx expected, P=%.2e %s)\n", 
                    sp_names[i], sp_names[j], overlap, 
                    test$fold_enrichment, test$p_value,
                    ifelse(test$p_value < 0.001, "***",
                           ifelse(test$p_value < 0.01, "**",
                                  ifelse(test$p_value < 0.05, "*", "ns")))))
      }
    }
  }
}

# Calculate three-way overlap for queens
if (length(queen_sets) == 3) {
  queen_all_three <- Reduce(intersect, queen_sets)
  cat(sprintf("  All three species: %d orthogroups\n", length(queen_all_three)))
}

cat("\nWorker-upregulated orthogroups:\n")
for (name in names(worker_sets)) {
  cat(sprintf("  %s: %d orthogroups\n", name, length(worker_sets[[name]])))
}

# Calculate pairwise overlaps for workers
if (length(worker_sets) >= 2) {
  sp_names <- names(worker_sets)
  for (i in 1:(length(sp_names)-1)) {
    for (j in (i+1):length(sp_names)) {
      overlap <- length(intersect(worker_sets[[sp_names[i]]], worker_sets[[sp_names[j]]]))
      test_name <- paste(species[i], "vs", species[j], "WORKER")
      if (test_name %in% names(overlap_tests)) {
        test <- overlap_tests[[test_name]]
        cat(sprintf("  %s & %s: %d orthogroups (%.2fx expected, P=%.2e %s)\n", 
                    sp_names[i], sp_names[j], overlap,
                    test$fold_enrichment, test$p_value,
                    ifelse(test$p_value < 0.001, "***",
                           ifelse(test$p_value < 0.01, "**",
                                  ifelse(test$p_value < 0.05, "*", "ns")))))
      }
    }
  }
}

# Calculate three-way overlap for workers
if (length(worker_sets) == 3) {
  worker_all_three <- Reduce(intersect, worker_sets)
  cat(sprintf("  All three species: %d orthogroups\n", length(worker_all_three)))
}

# ============================================================================
# 8. EXPORT RESULTS
# ============================================================================

# Save orthogroup assignments
write.csv(pattern_summary, "orthogroup_patterns.csv", row.names = FALSE)

# Save statistical tests (NEW)
write.csv(overlap_test_df, "overlap_significance_tests.csv", row.names = FALSE)

# Save conserved orthogroups
write.table(conserved_queen, "conserved_queen_orthogroups.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(conserved_worker, "conserved_worker_orthogroups.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(mixed_regulation, "mixed_regulation_orthogroups.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save intersection data
write.csv(all_int_df, "intersection_counts.csv", row.names = FALSE)

cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("\nGenerated files:\n")
cat("  FIGURES:\n")
cat("    - figure1_intersections.pdf/png (with significance stars)\n")
cat("    - figure2_heatmap.pdf\n")
cat("    - figure3_patterns.pdf/png\n")
cat("    - figure4_within_species.pdf/png\n")
cat("    - figure5_enrichment.pdf/png (NEW - fold enrichment analysis)\n")
cat("    - figure6_volcano.pdf/png (NEW - enrichment vs significance)\n")
cat("    - figure7_venn_queen.pdf/png (Venn diagram for queen orthogroups)\n")
cat("    - figure7_venn_worker.pdf/png (Venn diagram for worker orthogroups)\n")
cat("  DATA:\n")
cat("    - orthogroup_patterns.csv\n")
cat("    - overlap_significance_tests.csv (NEW - hypergeometric tests)\n")
cat("    - conserved_queen_orthogroups.txt\n")
cat("    - conserved_worker_orthogroups.txt\n")
cat("    - mixed_regulation_orthogroups.txt\n")
cat("    - intersection_counts.csv\n")
cat("\nSignificance codes: *** P<0.001, ** P<0.01, * P<0.05, ns P>=0.05\n")
cat("\nStatistical testing: Hypergeometric tests assess whether overlaps\n")
cat("are greater than expected by chance given the universe of orthogroups.\n")
