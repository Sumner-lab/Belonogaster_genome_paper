suppressMessages(library(DESeq2))
rc <- ""   # folder holding the count tables
spec <- list(BJ=c("Bjun_counts_table.csv",14,11), LF=c("LF_counts_table.csv",19,29), PD=c("Pdom_counts_table.csv",26,12))
pv <- function(m, sc=FALSE){ m <- m[matrixStats::rowVars(m)>0,]; p <- prcomp(t(m), scale.=sc); round(100*p$sdev^2/sum(p$sdev^2),2)[1:2] }
for (sp in names(spec)) { s <- spec[[sp]]
  cts <- as.matrix(read.csv(paste0(rc,s[1]), row.names=1, check.names=FALSE))
  dds <- DESeqDataSetFromMatrix(cts, data.frame(caste=factor(rep(c("rep","non"),c(as.integer(s[2]),as.integer(s[3])))),row.names=colnames(cts)), ~caste)
  dds <- estimateSizeFactors(dds); n <- counts(dds, normalized=TRUE)
  vb <- assay(vst(dds, blind=TRUE))
  cat(sp, "log2norm scaled", pv(log2(n+1),TRUE), "| vst scaled", pv(vb,TRUE), "| log10norm", pv(log10(n+1)), "| raw counts", pv(cts), "| norm counts", pv(n), "| CPM log2", pv(log2(t(t(cts)/colSums(cts))*1e6+1)), "\n") }
