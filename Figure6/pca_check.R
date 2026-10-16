suppressMessages(library(DESeq2))
rc <- ""   # folder holding the count tables
spec <- list(BJ=c("Bjun_counts_table.csv",14,11), LF=c("LF_counts_table.csv",19,29), PD=c("Pdom_counts_table.csv",26,12))
pv <- function(m, ntop) { v <- matrixStats::rowVars(m); sel <- order(v, decreasing=TRUE)[seq_len(min(ntop, length(v)))]
  p <- prcomp(t(m[sel,])); round(100*p$sdev^2/sum(p$sdev^2),2)[1:2] }
for (sp in names(spec)) { s <- spec[[sp]]
  cts <- as.matrix(read.csv(paste0(rc,s[1]), row.names=1, check.names=FALSE))
  caste <- factor(c(rep("rep",as.integer(s[2])),rep("non",as.integer(s[3]))))
  dds <- DESeqDataSetFromMatrix(cts, data.frame(caste=caste,row.names=colnames(cts)), ~caste)
  vb <- assay(vst(dds, blind=TRUE)); vf <- assay(vst(dds, blind=FALSE))
  rl <- assay(rlog(dds, blind=TRUE))
  dds <- estimateSizeFactors(dds); ln <- log2(counts(dds, normalized=TRUE)+1)
  cat(sp, "vst blind top500", pv(vb,500), "| vst all", pv(vb,1e6), "| vst nonblind top500", pv(vf,500),
      "| rlog top500", pv(rl,500), "| rlog all", pv(rl,1e6), "| log2norm top500", pv(ln,500), "| log2norm all", pv(ln,1e6), "\n") }
