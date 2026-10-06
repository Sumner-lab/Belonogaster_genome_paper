## ============================================================
## GO enrichment of caste-biased DEGs — B. juncea, L. flavolineata, P. dominula
##
## Run from this folder:  Rscript 01_run_topGO.R
##
## - DEGs: DESeq2 padj < 0.05 (BH), split by sign of log2FoldChange.
##   In all three DESeq2 tables log2FoldChange > 0 = higher in reproductives (queens),
##   log2FoldChange < 0 = higher in non-reproductives (workers).
## - Annotation: gene-level GO files from eggNOG-mapper 2.1.13, run by excon v2.5.0 on the
##   final annotations (*_final_annot.gff) with --tax_scope 50557 (Insecta) --score 60
##   --pident 40 --query_cover 40 --subject_cover 40  (go_files_excon_v2.5.0/).
## - Test: topGO "classic" algorithm + Fisher's exact test, BP / MF / CC separately.
##   Same test as the original ChopGO runs, but keeps ALL tested terms
##   (ChopGO kept only the top 50 per ontology and adjusted p-values over those 50).
## - Universe (background): expressed genes = genes with >= 1 GO annotation that are in
##   the DESeq2 table, i.e. that passed the low-expression filter used before DESeq2.
##   Robustness check: all genes with >= 1 GO annotation in that species ("global").
## - FDR: Benjamini-Hochberg across all tested terms, within each ontology.
## ============================================================

suppressWarnings(suppressMessages({
  library(topGO)
  library(writexl)
}))

PADJ_CUTOFF  <- 0.05   # DEG threshold (DESeq2 padj)
TABLE_P      <- 0.05   # terms with nominal P below this go into Supplementary Table 3
ONTOLOGIES   <- c("BP", "MF", "CC")

SPECIES <- list(
  BJ = list(name  = "B. juncea",
            res   = "../../res_Bjun 1.tsv",
            go    = "go_files_excon_v2.5.0/Belonogaster_juncea.go.txt"),
  LF = list(name  = "L. flavolineata",
            res   = "../../res_LF.tsv",
            go    = "go_files_excon_v2.5.0/Liostenogaster_flavolineata.go.txt"),
  PD = list(name  = "P. dominula",
            res   = "../../res_Pdom.tsv",
            go    = "go_files_excon_v2.5.0/Polistes_dominula.go.txt")
)

dir.create("gene_lists", showWarnings = FALSE)
dir.create("results",    showWarnings = FALSE)

## ── topGO for one gene set against one universe ──────────────────────────────

run_topgo <- function(gene2GO, sel, universe) {
  g2g <- gene2GO[names(gene2GO) %in% universe]
  inGenes <- factor(as.integer(names(g2g) %in% sel), levels = c(0, 1))
  names(inGenes) <- names(g2g)

  do.call(rbind, lapply(ONTOLOGIES, function(ont) {
    GOdata <- new("topGOdata", ontology = ont, allGenes = inGenes,
                  annot = annFUN.gene2GO, gene2GO = g2g)
    res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    p   <- score(res)
    tab <- GenTable(GOdata, classicFisher = res, orderBy = "classicFisher",
                    ranksOf = "classicFisher", topNodes = length(p), numChar = 1000)
    expected <- tab$Annotated * numSigGenes(GOdata) / numGenes(GOdata)

    out <- data.frame(
      GO_ID           = tab$GO.ID,
      GO_term         = tab$Term,
      Ontology        = ont,
      Annotated       = tab$Annotated,
      Significant     = tab$Significant,
      Expected        = round(expected, 3),
      Fold_enrichment = round(tab$Significant / expected, 3),
      P_value         = unname(p[tab$GO.ID]),
      stringsAsFactors = FALSE
    )
    out$FDR_BH <- p.adjust(out$P_value, method = "BH")

    # DEG IDs behind each nominally significant term
    out$DEG_IDs <- ""
    keep <- which(out$P_value < TABLE_P)
    if (length(keep)) {
      sig <- sigGenes(GOdata)
      git <- genesInTerm(GOdata, out$GO_ID[keep])
      out$DEG_IDs[keep] <- vapply(out$GO_ID[keep], function(id)
        paste(sort(intersect(git[[id]], sig)), collapse = ","), "")
    }
    out
  }))
}

## ── Export gene-to-term membership (as topGO sees it) for 03_go_overlap.py ──
## One file per species: every GO term in the topGO graph for each ontology and
## the background genes annotated to it (incl. inherited annotations), plus the
## background gene list. Used for the permutation test of GO-term overlap.

export_membership <- function(gene2GO, universe, sp) {
  g2g <- gene2GO[names(gene2GO) %in% universe]
  dummy <- factor(c(1, rep(0, length(g2g) - 1)), levels = c(0, 1))
  names(dummy) <- names(g2g)
  con <- gzfile(file.path("results", "membership", paste0(sp, "_membership.tsv.gz")), "w")
  writeLines("GO_ID\tOntology\tgene", con)
  for (ont in ONTOLOGIES) {
    GOdata <- new("topGOdata", ontology = ont, allGenes = dummy,
                  annot = annFUN.gene2GO, gene2GO = g2g)
    git <- genesInTerm(GOdata)
    writeLines(unlist(lapply(names(git), function(id)
      paste(id, ont, git[[id]], sep = "\t")), use.names = FALSE), con)
  }
  close(con)
  writeLines(names(g2g), file.path("results", "membership", paste0(sp, "_background.txt")))
}

dir.create(file.path("results", "membership"), showWarnings = FALSE)

## ── Run all species × caste ──────────────────────────────────────────────────

summary_rows  <- list()
supp_sheets   <- list()
all_results   <- list()
manifest_rows <- list()

for (sp in names(SPECIES)) {
  info <- SPECIES[[sp]]

  res <- read.delim(info$res, stringsAsFactors = FALSE)
  degs <- list(
    queen  = res$gene[!is.na(res$padj) & res$padj < PADJ_CUTOFF & res$log2FoldChange > 0],
    worker = res$gene[!is.na(res$padj) & res$padj < PADJ_CUTOFF & res$log2FoldChange < 0]
  )

  go <- read.table(info$go, sep = "\t", quote = "", comment.char = "",
                   col.names = c("gene", "go"), stringsAsFactors = FALSE)
  gene2GO <- lapply(split(go$go, go$gene), unique)
  export_membership(gene2GO, universe = res$gene, sp = sp)

  for (caste in names(degs)) {
    set_id <- paste(sp, caste, sep = "_")
    sel <- degs[[caste]]
    writeLines(sel, file.path("gene_lists", paste0(set_id, "_biased.txt")))

    # Main: universe = expressed genes (in the DESeq2 table, after expression filtering)
    main <- run_topgo(gene2GO, sel, universe = res$gene)
    write.table(main, file.path("results", paste0(set_id, "_biased_all_terms.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    all_results[[set_id]] <- cbind(species = sp, caste = caste, main)

    # Robustness: universe = all GO-annotated genes in the genome
    robust <- run_topgo(gene2GO, sel, universe = names(gene2GO))
    top_main   <- head(main$GO_ID[order(main$P_value)][main$P_value[order(main$P_value)] < 0.01], 12)
    top_robust <- head(robust$GO_ID[order(robust$P_value)][robust$P_value[order(robust$P_value)] < 0.01], 12)

    summary_rows[[set_id]] <- data.frame(
      species = sp, caste = caste,
      DEGs = length(sel),
      DEGs_with_GO = sum(sel %in% names(gene2GO)),
      background_genes = sum(names(gene2GO) %in% res$gene),
      terms_tested = nrow(main),
      terms_P_lt_0.01 = sum(main$P_value < 0.01),
      terms_FDR_lt_0.05 = sum(main$FDR_BH < 0.05),
      global_background_genes = length(gene2GO),
      global_terms_P_lt_0.01 = sum(robust$P_value < 0.01),
      global_terms_FDR_lt_0.05 = sum(robust$FDR_BH < 0.05),
      global_shared_P_lt_0.01 = length(intersect(main$GO_ID[main$P_value < 0.01],
                                                 robust$GO_ID[robust$P_value < 0.01])),
      global_shared_FDR_lt_0.05 = length(intersect(main$GO_ID[main$FDR_BH < 0.05],
                                                   robust$GO_ID[robust$FDR_BH < 0.05])),
      global_top12_shared = length(intersect(top_main, top_robust))
    )

    tab <- main[main$P_value < TABLE_P, ]
    tab <- tab[order(tab$P_value), ]
    list_name <- paste(info$name, paste0(caste, "-biased"))
    supp_sheets[[list_name]] <- tab
    # paths relative to results/ (where the manifest is written)
    manifest_rows[[set_id]] <- data.frame(
      list = list_name, species = sp, caste = caste,
      results_file    = paste0(set_id, "_biased_all_terms.tsv"),
      degs_file       = file.path("..", "gene_lists", paste0(set_id, "_biased.txt")),
      membership_file = file.path("membership", paste0(sp, "_membership.tsv.gz")),
      background_file = file.path("membership", paste0(sp, "_background.txt"))
    )
    message(sprintf("%s: %d DEGs (%d with GO), %d terms P<0.01, %d terms FDR<0.05",
                    set_id, length(sel), sum(sel %in% names(gene2GO)),
                    sum(main$P_value < 0.01), sum(main$FDR_BH < 0.05)))
  }
}

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, "results/summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, manifest_rows), "results/overlap_manifest.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_df, row.names = FALSE)

## ── Cross-species overlap of nominally enriched terms (P < 0.01), per caste ──

combined <- do.call(rbind, all_results)
write.table(combined, "results/all_sets_all_terms.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nTerms with P < 0.01 shared across species (same caste):\n")
for (caste in c("queen", "worker")) {
  sets <- lapply(names(SPECIES), function(sp)
    combined$GO_ID[combined$species == sp & combined$caste == caste & combined$P_value < 0.01])
  names(sets) <- names(SPECIES)
  pairs <- combn(names(sets), 2)
  for (j in seq_len(ncol(pairs))) {
    a <- pairs[1, j]; b <- pairs[2, j]
    shared <- intersect(sets[[a]], sets[[b]])
    cat(sprintf("  %s-biased  %s & %s: %d shared%s\n", caste, a, b, length(shared),
                if (length(shared)) paste0(" (", paste(shared, collapse = ", "), ")") else ""))
  }
}

## ── Supplementary Table 3 ────────────────────────────────────────────────────

about <- data.frame(
  Column = c("GO_ID", "GO_term", "Ontology", "Annotated", "Significant", "Expected",
             "Fold_enrichment", "P_value", "FDR_BH", "DEG_IDs"),
  Description = c(
    "Gene Ontology identifier",
    "GO term name",
    "BP = biological process, MF = molecular function, CC = cellular component",
    "Number of genes in the background (expressed genes, i.e. retained after low-expression filtering for DESeq2, with an eggNOG-derived GO annotation) annotated to the term, including annotations inherited from child terms",
    "Number of caste-biased DEGs annotated to the term",
    "Number of caste-biased DEGs expected to be annotated to the term by chance",
    "Significant / Expected",
    "Nominal P-value, one-sided Fisher's exact test (topGO 'classic' algorithm)",
    "Benjamini-Hochberg adjusted P-value across all GO terms tested within that ontology",
    "IDs of the caste-biased DEGs annotated to the term"
  )
)
about <- rbind(
  data.frame(Column = "Contents",
             Description = paste0("GO terms enriched (nominal P < ", TABLE_P,
                                  ") among queen-biased (higher in reproductives) and worker-biased ",
                                  "(higher in non-reproductives) DEGs (DESeq2 FDR < ", PADJ_CUTOFF,
                                  ") in each species. One sheet per species and caste.")),
  about
)

write_xlsx(c(list(About = about), supp_sheets), "Supplementary_Table_3_GO_enrichment.xlsx")
message("Saved: Supplementary_Table_3_GO_enrichment.xlsx")
