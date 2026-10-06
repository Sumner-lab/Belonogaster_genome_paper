## ============================================================
## GO enrichment dot plot (Figure 7) — adapted from ../R_plot/enrich.7.R
##
## Run from this folder after 01_run_topGO.R:  Rscript 02_plot_GO_figure.R
##
## Rows = species, columns = caste bias. Dot size = fold enrichment (capped at 20),
## fill = ontology, x = -log10(nominal Fisher P). Terms that also pass
## BH FDR < 0.05 (across all terms tested in that ontology) get a black outline.
## ============================================================

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(stringr)
  library(ggh4x)
}))

## ── 0. Configuration ──────────────────────────────────────────────────────────

FILES <- list(
  list(file = "results/BJ_queen_biased_all_terms.tsv",  species = "Belonogaster",   caste = "queen"),
  list(file = "results/BJ_worker_biased_all_terms.tsv", species = "Belonogaster",   caste = "worker"),
  list(file = "results/PD_queen_biased_all_terms.tsv",  species = "Polistes",       caste = "queen"),
  list(file = "results/PD_worker_biased_all_terms.tsv", species = "Polistes",       caste = "worker"),
  list(file = "results/LF_queen_biased_all_terms.tsv",  species = "Liostenogaster", caste = "queen"),
  list(file = "results/LF_worker_biased_all_terms.tsv", species = "Liostenogaster", caste = "worker")
)

CASTE_LABELS <- c(queen = "Queen-biased", worker = "Worker-biased")
PVAL_CUTOFF  <- 0.01   # nominal P threshold for a term to be shown
FDR_CUTOFF   <- 0.05   # terms below this get a black outline
TOP_N        <- 12     # top N terms per species × caste panel
ONTOLOGIES   <- c("BP", "MF", "CC")

## ── 1. Read and combine data ──────────────────────────────────────────────────

read_go <- function(info) {
  df <- read.delim(info$file, quote = "", stringsAsFactors = FALSE)
  df$species <- info$species
  df$caste   <- CASTE_LABELS[[info$caste]]
  df
}

all_go <- bind_rows(lapply(FILES, read_go))

all_go$ShortTerm <- str_wrap(all_go$GO_term, width = 30)
all_go$species   <- factor(all_go$species,
                           levels = c("Belonogaster", "Polistes", "Liostenogaster"))
all_go$caste     <- factor(all_go$caste, levels = unname(CASTE_LABELS))
all_go$Ontology  <- factor(all_go$Ontology, levels = ONTOLOGIES)

## ── 2. Filter and select top terms ────────────────────────────────────────────

top_terms <- all_go %>%
  filter(P_value < PVAL_CUTOFF) %>%
  mutate(neg_log10_p = -log10(P_value),
         fdr_sig = factor(if_else(FDR_BH < FDR_CUTOFF,
                                  paste0("FDR < ", FDR_CUTOFF), "Nominal P only"),
                          levels = c(paste0("FDR < ", FDR_CUTOFF), "Nominal P only"))) %>%
  group_by(species, caste) %>%
  slice_min(order_by = P_value, n = TOP_N, with_ties = FALSE) %>%
  mutate(
    # Unique y key per panel; label function strips the prefix
    panel_id = paste(species, caste, sep = "|||"),
    label_id = paste(panel_id, ShortTerm, sep = "~~~"),
    label_id = fct_reorder(label_id, desc(P_value))
  ) %>%
  ungroup()

## ── 3. Dot plot ───────────────────────────────────────────────────────────────

onto_colours  <- c(BP = "#4E79A7", MF = "#F28E2B", CC = "#59A14F")
# Figure 6 convention (as in 03_go_overlap.py): blue = reproductive, red = non-reproductive
caste_colours <- setNames(c("#2f6db5", "#c8423b"), unname(CASTE_LABELS))
fdr_outline   <- setNames(c("black", "white"), levels(top_terms$fdr_sig))

p_dot <- ggplot(top_terms,
                aes(x      = neg_log10_p,
                    y      = label_id,
                    fill   = Ontology,
                    colour = fdr_sig,
                    size   = pmin(Fold_enrichment, 20))) +
  geom_point(shape = 21, stroke = 0.6, alpha = 0.9) +
  ggh4x::facet_grid2(
    species ~ caste,
    scales      = "free_y",
    space       = "free_y",
    independent = "y",
    switch      = "y"
  ) +
  scale_y_discrete(labels = function(x) sub("^.+?~~~", "", x)) +
  scale_fill_manual(values = onto_colours, name = "Ontology", drop = FALSE) +
  scale_colour_manual(values = fdr_outline, name = NULL, drop = FALSE) +
  scale_size_continuous(
    range  = c(2, 9),
    name   = "Fold enrichment\n(capped at 20)",
    breaks = c(2, 5, 10, 20)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5, shape = 21, colour = "white", stroke = 0.4),
      order = 1
    ),
    size = guide_legend(
      override.aes = list(fill = "grey40", shape = 21, colour = "white", stroke = 0.4),
      order = 2
    ),
    colour = guide_legend(
      override.aes = list(fill = "grey70", shape = 21, size = 5, stroke = 0.8),
      order = 3
    )
  ) +
  labs(
    x = expression(-log[10](italic(p)~value)),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text.x       = element_text(size = 11, face = "bold", colour = "white"),
    strip.background.x = element_rect(fill = "grey35"),   # placeholder; recoloured below
    strip.text.y.left  = element_text(size = 9, angle = 0, hjust = 1, face = "italic"),
    strip.background.y = element_rect(fill = "grey85", colour = NA),
    strip.placement    = "outside",
    axis.text.y        = element_text(size = 7.5, lineheight = 0.85),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "right"
  )

# Recolour the caste column strip backgrounds
colour_caste_strips <- function(p) {
  gt <- ggplot_gtable(ggplot_build(p))
  strip_t <- which(grepl("strip-t", gt$layout$name))
  caste_levels <- levels(top_terms$caste)
  for (i in seq_along(strip_t)) {
    d <- caste_levels[((i - 1) %% length(caste_levels)) + 1]
    strip_grob <- gt$grobs[[strip_t[i]]]
    tryCatch({
      strip_grob$grobs[[1]]$children[[1]]$gp$fill <- caste_colours[d]
      gt$grobs[[strip_t[i]]] <- strip_grob
    }, error = function(e) NULL)
  }
  gt
}

p_dot_gt <- colour_caste_strips(p_dot)

## ── 4. Save ───────────────────────────────────────────────────────────────────

pdf("Figure7_GO_dotplot_caste.pdf", width = 14, height = 16)
grid::grid.draw(p_dot_gt)
dev.off()

png("Figure7_GO_dotplot_caste.png", width = 14, height = 16, units = "in", res = 150)
grid::grid.draw(p_dot_gt)
dev.off()

message("Saved: Figure7_GO_dotplot_caste.pdf / .png")
