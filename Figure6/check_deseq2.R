suppressMessages(library(DESeq2))
rc <- ""   # folder holding the count tables
spec <- list(
  BJ = list(f="Bjun_counts_table.csv", res="res_Bjun 1.tsv", nrep=14, nnon=11),
  LF = list(f="LF_counts_table.csv",   res="res_LF.tsv",     nrep=19, nnon=29),
  PD = list(f="Pdom_counts_table.csv", res="res_Pdom.tsv",   nrep=26, nnon=12))
for (sp in names(spec)) {
  s <- spec[[sp]]
  cts <- as.matrix(read.csv(paste0(rc, s$f), row.names = 1, check.names = FALSE))
  caste <- factor(c(rep("rep", s$nrep), rep("non", s$nnon)), levels = c("non", "rep"))
  dds <- DESeqDataSetFromMatrix(cts, data.frame(caste = caste, row.names = colnames(cts)), ~ caste)
  dds <- DESeq(dds, fitType = "local", quiet = TRUE)
  r <- as.data.frame(results(dds)); r$gene <- rownames(r)
  old <- read.delim(s$res)
  m <- merge(old, r, by = "gene", suffixes = c("_old", "_new"))
  cat(sprintf("%s: %d genes in counts, %d in old table; LFC cor = %.6f; max |dLFC| = %.2g; DEGs old %d new %d\n",
      sp, nrow(cts), nrow(old), cor(m$log2FoldChange_old, m$log2FoldChange_new, use = "complete"),
      max(abs(m$log2FoldChange_old - m$log2FoldChange_new), na.rm = TRUE),
      sum(old$padj < 0.05, na.rm = TRUE), sum(r$padj < 0.05, na.rm = TRUE)))
}
