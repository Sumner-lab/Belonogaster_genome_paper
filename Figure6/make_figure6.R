#!/usr/bin/env Rscript
# Figure 6: differential expression between castes in B. juncea, L. flavolineata and P. dominula.
#  top    : PCA of variance-stabilised counts (DESeq2 vst, blind = TRUE; 500 most variable genes, as plotPCA)
#  middle : volcano plots from DESeq2 (~ caste, fitType = "local"; reproductives vs non-reproductives)
#  bottom : pairwise overlaps of caste-biased orthogroups (excon run 13, OrthoFinder v3.1.4), with
#           one-sided hypergeometric P (universe = orthogroups with a tested gene in both species)
# The DESeq2 fits reproduce res_Bjun 1.tsv, res_LF.tsv and res_Pdom.tsv exactly (check_deseq2.R).
# Count columns are ordered reproductives first, then non-reproductives (Supp Table 2: Samples).
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(patchwork) })

rc <- ""   # folder holding the count tables
spec <- list(
  BJ = list(name = "B. juncea",       f = "Bjun_counts_table.csv", nrep = 14, nnon = 11),
  LF = list(name = "L. flavolineata", f = "LF_counts_table.csv",   nrep = 19, nnon = 29),
  PD = list(name = "P. dominula",     f = "Pdom_counts_table.csv", nrep = 26, nnon = 12))
col_status <- c(Nonreproductive = "#E41A1C", Reproductive = "#377EB8")
shp_status <- c(Nonreproductive = 16, Reproductive = 17)
base <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

pca_plots <- list(); volc_plots <- list(); deg_tab <- list()
for (sp in names(spec)) {
  s <- spec[[sp]]
  cts <- as.matrix(read.csv(paste0(rc, s$f), row.names = 1, check.names = FALSE))
  status <- factor(rep(c("Reproductive", "Nonreproductive"), c(s$nrep, s$nnon)),
                   levels = c("Nonreproductive", "Reproductive"))
  dds <- DESeqDataSetFromMatrix(cts, data.frame(caste = status, row.names = colnames(cts)), ~ caste)
  dds <- DESeq(dds, fitType = "local", quiet = TRUE)
  res <- as.data.frame(results(dds))          # log2FC = reproductive / non-reproductive
  deg_tab[[sp]] <- c(sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
                     sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE))

  pd <- plotPCA(vst(dds, blind = TRUE), intgroup = "caste", returnData = TRUE, ntop = 500)
  pvar <- round(100 * attr(pd, "percentVar"), 1)
  pca_plots[[sp]] <- ggplot(pd, aes(PC1, PC2, colour = caste, shape = caste)) +
    geom_point(size = 2.4, alpha = 0.9) +
    scale_colour_manual(values = col_status, name = "Reproductive status") +
    scale_shape_manual(values = shp_status, name = "Reproductive status") +
    labs(x = sprintf("PC1 (%.1f%%)", pvar[1]), y = sprintf("PC2 (%.1f%%)", pvar[2]), title = s$name) +
    base + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 12))

  res <- res[!is.na(res$padj), ]
  res$sig <- res$padj < 0.05
  volc_plots[[sp]] <- ggplot(res, aes(log2FoldChange, -log10(padj))) +
    geom_point(data = res[!res$sig, ], colour = "grey25", size = 0.8, alpha = 0.6) +
    geom_point(data = res[res$sig, ], colour = "#E41A1C", size = 0.8, alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
    labs(x = expression(log[2]~"fold change (rep / nonrep)"),
         y = expression(-log[10]~"FDR-adjusted"~italic(P))) + base
}
print(deg_tab)

# ---- bottom row: pairwise orthogroup overlaps (from overlap_test.py output) ----
ov <- read.delim("../DEG_overlap/overlap_tests_excon_run13_OFv3.tsv")
ov <- ov[ov$universe == "B_tested_both", ]
lab <- c(Queen = "rep", Worker = "nonrep")
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
getov <- function(a, ca, b, cb) {
  r <- ov[(ov$set1 == paste(a, ca) & ov$set2 == paste(b, cb)) | (ov$set1 == paste(b, cb) & ov$set2 == paste(a, ca)), ]
  c(n = r$overlap, p = r$P)
}
upset_pair <- function(a, b) {
  combos <- list(c("Worker", "Worker"), c("Queen", "Queen"), c("Worker", "Queen"), c("Queen", "Worker"))
  rows <- c(paste(a, "nonrep"), paste(b, "rep"), paste(b, "nonrep"), paste(a, "rep"))
  d <- do.call(rbind, lapply(combos, function(cc) {
    o <- getov(a, cc[1], b, cc[2])
    data.frame(n = o["n"], s = if (cc[1] == cc[2]) star(o["p"]) else "",
               r1 = paste(a, lab[cc[1]]), r2 = paste(b, lab[cc[2]]))
  }))
  d <- d[order(-d$n), ]; d$x <- seq_len(nrow(d))
  top <- max(d$n)
  bars <- ggplot(d, aes(x, n)) + geom_col(fill = "grey20", width = 0.6) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.2) +
    geom_text(aes(y = n + top * 0.13, label = s), vjust = -0.2, size = 3.8, fontface = "bold") +
    scale_y_continuous(limits = c(0, top * 1.32), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0.4, 4.6), expand = c(0, 0)) +
    labs(x = NULL, y = "Shared orthogroups",
         title = paste(spec[[a]]$name, "vs", spec[[b]]$name)) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 11))
  grid <- expand.grid(x = d$x, row = rows)
  grid$on <- mapply(function(x, r) r %in% c(d$r1[d$x == x], d$r2[d$x == x]), grid$x, as.character(grid$row))
  seg <- data.frame(x = d$x, y1 = match(d$r1, rev(rows)), y2 = match(d$r2, rev(rows)))
  grid$y <- match(as.character(grid$row), rev(rows))
  dots <- ggplot() +
    annotate("rect", xmin = 0.45, xmax = 4.55, ymin = c(0.5, 2.5), ymax = c(1.5, 3.5), fill = "grey92") +
    geom_point(data = grid[!grid$on, ], aes(x, y), colour = "grey80", size = 1.6) +
    geom_segment(data = seg, aes(x = x, xend = x, y = y1, yend = y2), linewidth = 0.7) +
    geom_point(data = grid[grid$on, ], aes(x, y), size = 2.2) +
    scale_y_continuous(breaks = 1:4, labels = rev(rows), expand = c(0, 0.5)) +
    scale_x_continuous(limits = c(0.4, 4.6), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) + theme_void(base_size = 10) +
    theme(axis.text.y = element_text(size = 8.5, hjust = 1, margin = margin(r = 3)))
  bars / dots + plot_layout(heights = c(3, 1.3))
}
bottom <- wrap_plots(upset_pair("BJ", "LF"), upset_pair("BJ", "PD"), upset_pair("LF", "PD"), nrow = 1)

top <- wrap_plots(pca_plots, nrow = 1, guides = "collect") & theme(legend.position = "bottom")
mid <- wrap_plots(volc_plots, nrow = 1)
fig <- wrap_elements(top) / wrap_elements(mid) / wrap_elements(bottom) +
  plot_layout(heights = c(1.15, 1, 1.25)) +
  plot_annotation(tag_levels = list(c("a", "b", "c"))) &
  theme(plot.tag = element_text(face = "bold", size = 14))
qpdf <- function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
ggsave("Figure6.pdf", fig, width = 10, height = 11, device = qpdf)
ggsave("Figure6.png", fig, width = 10, height = 11, dpi = 300, bg = "white")
ggsave("Figure6.tif", fig, width = 10, height = 11, dpi = 300, bg = "white", compression = "lzw")
cat("Wrote Figure6.pdf/.png/.tif\n")
