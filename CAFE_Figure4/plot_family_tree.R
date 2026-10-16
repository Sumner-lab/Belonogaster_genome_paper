#!/usr/bin/env Rscript
###############################################################################
# plot_family_tree.R: copy number of one gene family on every branch of the species tree
# (Figure 4c). Same tree layout and tip order as plot_go_tree_heatmap.R --cladogram.
#
# Usage:
#   LC_ALL=en_GB.UTF-8 Rscript plot_family_tree.R species_tree.nwk cafe_numbered_tree.nwk \
#       CAFE_dir FAMILY_ID out_prefix [grey_branches]
#   CAFE_dir holds Gamma_count.tab, Gamma_change.tab and Gamma_branch_probabilities.tab.
#   grey_branches: comma-separated tips / CAFE node ids drawn grey (uninterpretable root branches).
###############################################################################
suppressPackageStartupMessages({ library(ape); library(ggtree); library(ggplot2); library(dplyr) })

a <- commandArgs(trailingOnly = TRUE)
tree <- read.tree(a[1]); tc <- read.tree(a[2]); cdir <- a[3]; fam <- a[4]; out <- a[5]
tree$tip.label <- sub("\\.clean$", "", tree$tip.label); tc$tip.label <- sub("\\.clean$", "", tc$tip.label)
grey <- if (length(a) >= 6) trimws(strsplit(a[6], ",")[[1]]) else character(0)

row_of <- function(f) {
  x <- read.delim(file.path(cdir, f), check.names = FALSE, colClasses = "character")
  r <- unlist(x[x[[1]] == fam, -1])
  names(r) <- sub(">$", "", sub("^.*<", "", names(x)[-1]))   # "Species<3>" -> "3"
  r
}
cnt <- row_of("Gamma_count.tab"); chg <- row_of("Gamma_change.tab"); prb <- row_of("Gamma_branch_probabilities.tab")

# CAFE id of every tip and internal node of `tree`
x <- read.delim(file.path(cdir, "Gamma_count.tab"), check.names = FALSE, nrows = 1)
tipid <- setNames(sub(">$", "", sub("^.*<", "", names(x)[-1])), sub("<.*$", "", names(x)[-1]))
tipid <- tipid[nzchar(names(tipid))]
d <- as.data.frame(ggtree(tree, branch.length = "none")$data)
d$cafe <- NA_character_
d$cafe[d$isTip] <- tipid[d$label[d$isTip]]
ntt <- length(tc$tip.label)
for (i in seq_len(tc$Nnode)) {
  ts <- extract.clade(tc, ntt + i)$tip.label
  n <- if (length(ts) == ntt) length(tree$tip.label) + 1L else getMRCA(tree, ts)
  d$cafe[d$node == n] <- tc$node.label[i]
}
d$parent_x <- d$x[match(d$parent, d$node)]
d$count <- cnt[d$cafe]; d$change <- as.integer(chg[d$cafe])
d$p <- suppressWarnings(as.numeric(prb[d$cafe]))
d$stars <- ifelse(is.na(d$p) | d$p >= 0.05, "", ifelse(d$p < 0.001, "***", ifelse(d$p < 0.01, "**", "*")))
d$col <- ifelse(d$stars == "", "black", ifelse(d$change > 0, "#1b7837", "#b2182b"))
d$col[d$cafe %in% grey | d$label %in% grey] <- "grey60"
d$lab <- paste0(d$stars, d$count)

seg_h <- data.frame(x = d$parent_x, xend = d$x, y = d$y)
seg_v <- do.call(rbind, lapply(unique(d$parent), function(p) {
  k <- d[d$parent == p & d$node != p, ]
  if (!nrow(k)) return(NULL)
  data.frame(x = d$x[d$node == p][1], ymin = min(k$y), ymax = max(k$y))
}))
xmax <- max(d$x)

g <- ggplot() +
  geom_segment(data = seg_h, aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.4) +
  geom_segment(data = seg_v, aes(x = x, xend = x, y = ymin, yend = ymax), linewidth = 0.4) +
  geom_text(data = d[d$isTip, ], aes(x = x, y = y, label = gsub("_", " ", label)),
            hjust = 0, nudge_x = xmax * 0.03, size = 3.6, fontface = "italic") +
  geom_text(data = d, aes(x = x - xmax * 0.02, y = y + 0.28, label = lab), colour = d$col,
            hjust = 1, size = 3.2, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.62))) +
  scale_y_continuous(limits = c(0.5, length(tree$tip.label) + 0.6), expand = c(0, 0)) +
  labs(caption = "*P < 0.05, **P < 0.01, ***P < 0.001; green = significant gain, red = significant loss on the branch") +
  theme_void() + theme(plot.caption = element_text(size = 8, hjust = 0), plot.margin = margin(4, 4, 4, 4))

qpdf <- function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
ggsave(paste0(out, ".pdf"), g, width = 5.5, height = 5, device = qpdf)
ggsave(paste0(out, ".png"), g, width = 5.5, height = 5, dpi = 600, bg = "white")
print(d[order(d$y), c("label", "cafe", "count", "change", "p", "stars")], row.names = FALSE)
cat("Wrote", out, ".pdf/.png\n")
