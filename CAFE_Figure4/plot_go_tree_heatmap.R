#!/usr/bin/env Rscript

# NOTE ON VERSIONS -------------------------------------------------------------
# ggtree 3.12 + ggplot2 4.0 are incompatible (ggtree's drawing geoms call the
# removed is.waive()). This script reads ONLY the layout coordinates from
# ggtree(tree)$data and draws the tree with plain geom_segment().
# Full GO names come from GO.db if installed; CAFE expansion/contraction counts
# from CAFE_summary.txt are drawn on the tree (leaves by name, internal nodes by
# CAFE id -- verify the printed mapping table against cafe_node_label_guide.pdf).
# ------------------------------------------------------------------------------

library(ape)
library(ggtree)
library(ggplot2)
library(tidyverse)

if(!requireNamespace("patchwork", quietly = TRUE)){
  stop("This script needs the 'patchwork' package:  install.packages('patchwork')")
}
library(patchwork)

###############################################################################
# ARGUMENTS  (positional + --flags)
###############################################################################
#
# Usage:
#   Rscript plot_go_tree_heatmap.R tree.nwk go.tsv [prefix] [cafe_summary.txt] [flags]
#
# Flags (all optional; --key value  or  --key=value):
#   --top N          max GO terms to show (default 50; use 'all' for no cap)
#   --min-species N  keep terms significant in >= N species (default 1)
#   --pval P         significance threshold (default 0.05)
#   --mask           colour only significant cells (off by default)
#   --label-size N   GO name / GO id axis font size (default 6)
#   --tip-size N     species-name font size (default 3)
#   --ann-size N     +exp / -con node-label font size (default 2.3)
#   --tree-frac N    tree-panel width relative to --heat-frac (default 1)
#   --heat-frac N    heatmap-panel width (default 3)
#   --label-pad N    right-side space for tip labels, fraction of tree (default 0.55)
#   --cladogram      draw all branches the same length, so +/- labels on short branches don't overlap
#   --grey-branches  comma-separated tips / CAFE node ids whose +/- labels are drawn grey (e.g. uninterpretable root branches)
#   --caption TEXT   caption under the figure (default: explains the +/- labels; 'none' to drop)
#
# e.g.  Rscript plot_go_tree_heatmap.R tree.nwk go.tsv wasp_GO CAFE_summary.txt \
#         --top 10 --tip-size 4 --label-size 8 --tree-frac 2
###############################################################################

raw   <- commandArgs(trailingOnly = TRUE)
flags <- list()
pos   <- character(0)
i <- 1
while(i <= length(raw)){
  a <- raw[i]
  if(grepl("^--", a)){
    if(grepl("=", a)){
      key <- sub("=.*$", "", sub("^--", "", a))
      flags[[key]] <- sub("^[^=]*=", "", a)
      i <- i + 1
    } else if(i + 1 <= length(raw) && !grepl("^--", raw[i + 1])){
      flags[[sub("^--", "", a)]] <- raw[i + 1]; i <- i + 2
    } else {
      flags[[sub("^--", "", a)]] <- "TRUE"; i <- i + 1   # boolean switch
    }
  } else {
    pos <- c(pos, a); i <- i + 1
  }
}
getf <- function(name, default) if(!is.null(flags[[name]])) flags[[name]] else default

if(length(pos) < 2){
  stop("Usage: Rscript plot_go_tree_heatmap.R tree.nwk go.tsv [prefix] [cafe.txt] [--flags]\n",
       "See the header comment for the list of flags.")
}
tree_file <- pos[1]
go_file   <- pos[2]
prefix    <- ifelse(length(pos) >= 3, pos[3], "GO_tree")
cafe_file <- ifelse(length(pos) >= 4, pos[4], "CAFE_summary.txt")

###############################################################################
# TUNABLES  (flag value, else default)
###############################################################################

PVAL_CUTOFF <- as.numeric(getf("pval", "0.05"))
MIN_SPECIES <- as.integer(getf("min-species", "1"))
TOP_N       <- { v <- getf("top", "50")
                 if(tolower(v) %in% c("na", "all", "none", "0")) NA_integer_ else as.integer(v) }
MASK_NONSIG <- tolower(getf("mask", "false")) %in% c("true", "1", "yes")
CLADOGRAM   <- tolower(getf("cladogram", "false")) %in% c("true", "1", "yes")
GREY_BR     <- trimws(strsplit(getf("grey-branches", ""), ",")[[1]])

LABEL_SIZE  <- as.numeric(getf("label-size", "6"))    # GO axis labels
TIP_SIZE    <- as.numeric(getf("tip-size",   "3"))    # species names
ANN_SIZE    <- as.numeric(getf("ann-size",   "2.3"))  # +exp / -con

TREE_FRAC   <- as.numeric(getf("tree-frac", "1"))
HEAT_FRAC   <- as.numeric(getf("heat-frac", "3"))
LABEL_PAD   <- as.numeric(getf("label-pad", "0.55"))

# --- internal-node mapping (edit these in-file if needed) --------------------
CAFE_TREE         <- getf("cafe-tree", "")    # newick whose internal labels ARE the CAFE numbers
FILL_LABEL        <- gsub("\\\\n", "\n", getf("fill-label", "-log10(p)"))
NODE_MAP_OVERRIDE <- c()   # e.g. c("18" = 19); names = CAFE id, values = ape node

CAPTION <- getf("caption", "Tree labels: green +n = significantly expanded gene families; red \u2212n = significantly contracted (CAFE, p \u2264 0.05)")
if(tolower(CAPTION) == "none") CAPTION <- NULL

EXP_COL <- "#1b7837"
CON_COL <- "#b2182b"

cat(sprintf("Settings: top=%s  min-species=%d  pval=%g  tip-size=%g  label-size=%g  tree:heat=%g:%g\n",
            ifelse(is.na(TOP_N), "all", TOP_N), MIN_SPECIES, PVAL_CUTOFF,
            TIP_SIZE, LABEL_SIZE, TREE_FRAC, HEAT_FRAC))

###############################################################################
# TREE
###############################################################################

cat("Reading tree\n")
tree <- read.tree(tree_file)
tree$tip.label <- sub("\\.clean$", "", tree$tip.label)
n_tips <- length(tree$tip.label)

###############################################################################
# GO TABLE
###############################################################################

cat("Reading GO table\n")
go <- read.delim(go_file, check.names = FALSE, stringsAsFactors = FALSE)
species_cols <- setdiff(names(go), c("GO_ID", "GO_term", "Count_significant"))
go <- go %>% distinct(GO_ID, .keep_all = TRUE)

if(requireNamespace("GO.db", quietly = TRUE)){
  lk <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      GO.db::GO.db, keys = unique(go$GO_ID), columns = "TERM", keytype = "GOID")),
    error = function(e) NULL
  )
  if(!is.null(lk)){
    full <- lk$TERM[match(go$GO_ID, lk$GOID)]
    go$GO_term <- ifelse(is.na(full) | full == "", go$GO_term, full)
    cat(sprintf("Recovered %d full GO term name(s) from GO.db\n",
                sum(!is.na(full) & full != "")))
  }
} else {
  cat("\n[note] GO.db not installed; GO names may be truncated.\n",
      "       Install with: BiocManager::install('GO.db')\n", sep = "")
}

pmat <- sapply(go[species_cols], as.numeric)
go$n_sig <- rowSums(!is.na(pmat) & pmat < PVAL_CUTOFF)

n_before <- nrow(go)
go <- go %>% filter(n_sig >= MIN_SPECIES)
if(!is.na(TOP_N) && nrow(go) > TOP_N){
  go <- go %>% arrange(desc(n_sig)) %>% slice_head(n = TOP_N)
}
cat(sprintf("\nKept %d of %d GO terms (MIN_SPECIES=%d, TOP_N=%s)\n",
            nrow(go), n_before, MIN_SPECIES,
            ifelse(is.na(TOP_N), "all", as.character(TOP_N))))
if(nrow(go) == 0) stop("No GO terms passed the filters.")

id_name <- setNames(go$GO_term, go$GO_ID)

score_long <- go %>%
  select(GO_ID, all_of(species_cols)) %>%
  pivot_longer(-GO_ID, names_to = "species", values_to = "pvalue") %>%
  mutate(
    species = sub("\\.clean$", "", species),
    pvalue  = as.numeric(pvalue),
    sig     = !is.na(pvalue) & pvalue < PVAL_CUTOFF,
    score   = ifelse(is.na(pvalue), NA_real_, -log10(pvalue))
  )
if(MASK_NONSIG) score_long <- score_long %>% mutate(score = ifelse(sig, score, NA_real_))

###############################################################################
# TREE COORDINATES -> manual segments
###############################################################################

d <- as.data.frame(ggtree(tree, branch.length = if(CLADOGRAM) "none" else "branch.length")$data)
d$parent_x <- d$x[match(d$parent, d$node)]
d$midx     <- (d$parent_x + d$x) / 2

seg_h <- data.frame(x = d$parent_x, xend = d$x, y = d$y, yend = d$y)
seg_v <- do.call(rbind, lapply(unique(d$parent), function(p){
  kids <- d[d$parent == p & d$node != p, , drop = FALSE]
  if(nrow(kids) < 1) return(NULL)
  px <- d$x[d$node == p][1]
  data.frame(x = px, xend = px, y = min(kids$y), yend = max(kids$y))
}))
tips      <- d[d$isTip, c("x", "y", "label")]
tree_xmax <- max(d$x)

###############################################################################
# CAFE EXPANSION / CONTRACTION ANNOTATIONS
###############################################################################

ann <- NULL

if(file.exists(cafe_file)){
  cat("\nReading CAFE summary\n")
  cafe <- read.delim(cafe_file, comment.char = "#",
                     check.names = FALSE, stringsAsFactors = FALSE)
  names(cafe)[1] <- "id"
  cafe$id <- trimws(as.character(cafe$id))

  need <- c("Expansion_HOGs_significant", "Contraction_HOGs_significant")
  if(!all(need %in% names(cafe))){
    cat("[warn] CAFE file missing expected columns; skipping annotations.\n")
  } else {
    cafe$exp <- as.integer(cafe[["Expansion_HOGs_significant"]])
    cafe$con <- as.integer(cafe[["Contraction_HOGs_significant"]])

    tip_ann <- d %>%
      filter(isTip) %>%
      inner_join(cafe %>% select(id, exp, con), by = c("label" = "id"))

    cafe_int <- cafe %>%
      filter(grepl("^[0-9]+$", id)) %>%
      mutate(cafe_id = as.integer(id))

    desc_tips <- function(n){
      if(n <= n_tips) return(tree$tip.label[n])
      tryCatch(extract.clade(tree, n)$tip.label, error = function(e) NA_character_)
    }

    if(nzchar(CAFE_TREE) && file.exists(CAFE_TREE)){
      cat("Deriving node map from", CAFE_TREE, "by descendant tips\n")
      tc <- read.tree(CAFE_TREE)
      tc$tip.label <- sub("\\.clean$", "", tc$tip.label)
      ntt <- length(tc$tip.label)
      m <- sapply(seq_len(tc$Nnode), function(i){
        cid <- tc$node.label[i]
        if(is.na(cid) || cid == "") return(NA_integer_)
        ts <- extract.clade(tc, ntt + i)$tip.label
        if(length(ts) == ntt) n_tips + 1L else getMRCA(tree, ts)
      })
      names(m) <- tc$node.label[seq_len(tc$Nnode)]
      cafe_int$ape_node <- m[as.character(cafe_int$cafe_id)]
    } else {
      cafe_int$ape_node <- cafe_int$cafe_id
      if(length(NODE_MAP_OVERRIDE)){
        hit <- as.character(cafe_int$cafe_id) %in% names(NODE_MAP_OVERRIDE)
        cafe_int$ape_node[hit] <-
          as.integer(NODE_MAP_OVERRIDE[as.character(cafe_int$cafe_id)[hit]])
      }
    }

    node_ann <- d %>%
      filter(!isTip) %>%
      inner_join(cafe_int %>% select(cafe_id, ape_node, exp, con),
                 by = c("node" = "ape_node"))

    if(nrow(node_ann) > 0){
      verify <- node_ann %>%
        transmute(cafe_id, ape_node = node, n_label = label, exp, con,
                  tips = sapply(node, function(n) paste(desc_tips(n), collapse = ", ")))
      cat("\n=== VERIFY internal-node mapping against cafe_node_label_guide.pdf ===\n")
      print(verify, right = FALSE, width = 250)
      cat("If a row's 'tips' don't match the CAFE node, set NODE_MAP_OVERRIDE",
          "or point CAFE_TREE at the numbered newick.\n")
    }

    unmatched <- setdiff(cafe_int$cafe_id, node_ann$cafe_id)
    if(length(unmatched) > 0){
      cat("\n[warn] CAFE internal nodes with no tree match:",
          paste(unmatched, collapse = ", "), "\n")
    }

    ann <- bind_rows(
      tip_ann  %>% transmute(id = label, midx, y, exp, con),
      node_ann %>% transmute(id = as.character(cafe_id), midx, y, exp, con)
    )
    ann$grey <- ann$id %in% GREY_BR
  }
} else {
  cat("\n[note] CAFE summary '", cafe_file,
      "' not found; tree drawn without expansion/contraction labels.\n", sep = "")
}

###############################################################################
# TREE PANEL
###############################################################################

p_tree <- ggplot() +
  geom_segment(data = seg_h, aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.4) +
  geom_segment(data = seg_v, aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.4) +
  geom_text(data = tips, aes(x = x, y = y, label = gsub("_", " ", label)),
            hjust = 0, size = TIP_SIZE, nudge_x = tree_xmax * 0.02, fontface = "italic")

if(!is.null(ann) && nrow(ann) > 0){
  p_tree <- p_tree +
    geom_text(data = ann, aes(x = midx, y = y + 0.24, label = paste0("+", exp)),
              colour = ifelse(ann$grey, "grey65", EXP_COL), size = ANN_SIZE, fontface = "bold") +
    geom_text(data = ann, aes(x = midx, y = y - 0.24, label = paste0("\u2212", con)),
              colour = ifelse(ann$grey, "grey65", CON_COL), size = ANN_SIZE, fontface = "bold")
}

p_tree <- p_tree +
  scale_y_continuous(limits = c(0.5, n_tips + 0.5), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, LABEL_PAD))) +
  theme_void() +
  theme(plot.margin = margin(6, 2, 6, 6))

tip_y <- tips %>% select(label, y)

###############################################################################
# HEATMAP PANEL
###############################################################################

make_heatmap <- function(ordered_ids){
  hm <- score_long %>%
    mutate(GO_ID = factor(GO_ID, levels = ordered_ids)) %>%
    left_join(tip_y, by = c("species" = "label")) %>%
    filter(!is.na(y))

  ggplot(hm, aes(x = GO_ID, y = y, fill = score)) +
    geom_tile(colour = "grey92", linewidth = 0.15) +
    scale_fill_gradient(low = "white", high = "firebrick",
                        na.value = "grey95", name = FILL_LABEL) +
    scale_y_continuous(limits = c(0.5, n_tips + 0.5), expand = c(0, 0)) +
    scale_x_discrete(limits = ordered_ids,
                     labels = unname(id_name[ordered_ids]),
                     sec.axis = dup_axis(labels = ordered_ids)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, size = LABEL_SIZE),
      axis.text.x.top    = element_text(angle = 90, hjust = 0, vjust = 0.5, size = LABEL_SIZE),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      plot.margin  = margin(6, 6, 6, 2)
    )
}

make_plot <- function(ordered_ids, file){
  combined <- p_tree + make_heatmap(ordered_ids) +
    plot_layout(widths = c(TREE_FRAC, HEAT_FRAC)) +
    plot_annotation(caption = CAPTION)
  # Quartz PDF (macOS) so the minus sign (U+2212) and <= render; run with a UTF-8 locale,
  # e.g. LC_ALL=en_GB.UTF-8 Rscript ... (cairo_pdf writes nothing on this machine). PNG copy for previews.
  qpdf <- function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
  ggsave(file, combined, width = 16, height = 14, limitsize = FALSE, device = qpdf)
  ggsave(sub("\\.pdf$", ".png", file), combined, width = 16, height = 14, dpi = 300, limitsize = FALSE)
  cat("Wrote", file, "and .png\n")
}

###############################################################################
# OUTPUT
###############################################################################

ord_sig <- go %>% arrange(desc(n_sig)) %>% pull(GO_ID)
make_plot(ord_sig, paste0(prefix, "_tree_heatmap.pdf"))

mat <- score_long %>%
  select(species, GO_ID, score) %>%
  pivot_wider(names_from = GO_ID, values_from = score) %>%
  as.data.frame()
rownames(mat) <- mat$species
mat$species   <- NULL

if(ncol(mat) > 2){
  cat("Building clustered version\n")
  m0 <- mat; m0[is.na(m0)] <- 0
  ord <- hclust(dist(t(m0)))$order
  make_plot(colnames(mat)[ord], paste0(prefix, "_tree_heatmap_clustered.pdf"))
} else {
  cat("Too few GO terms to cluster; skipping clustered plot\n")
}

cat("\nDone\n")
