# Figure 6: differential expression between castes

```
LC_ALL=en_GB.UTF-8 Rscript make_figure6.R
```
R 4.4 with DESeq2 1.44.0, ggplot2 and patchwork. The PDF uses the Quartz device (macOS); on other systems change `qpdf` to `cairo_pdf`.

## Inputs (put them in this folder)

- `Bjun_counts_table.csv`, `LF_counts_table.csv`, `Pdom_counts_table.csv`: brain RNA-seq gene counts, already expression-filtered (10,243, 10,934 and 9,254 genes). Columns are ordered reproductives first (14, 19, 26), then non-reproductives (11, 29, 12).
- `../DEG_overlap/overlap_tests_excon_run13_OFv3.tsv`: output of `DEG_overlap/overlap_test.py`, for panel c.

## Scripts

| Script | What it does |
|---|---|
| `make_figure6.R` | DESeq2 (`~ caste`, `fitType = "local"`); PCA, volcano plots and orthogroup overlaps; writes `Figure6.pdf/.png/.tif` |
| `check_deseq2.R` | Checks that this fit reproduces the published result tables (`res_Bjun 1.tsv`, `res_LF.tsv`, `res_Pdom.tsv`, Supplementary Table 2): fold changes agree to 1e-8, and the DEG counts are 573, 583 and 2,416 |
| `pca_check.R`, `pca_check2.R` | Compare PCA settings (vst / rlog, blind or not, number of genes) |

## Panels

- **a. PCA:** `vst(blind = TRUE)` and `plotPCA` defaults (500 most variable genes). PC1 / PC2: *B. juncea* 32.9 / 13.7%, *L. flavolineata* 6.2 / 5.3%, *P. dominula* 23.4 / 10.0%.
- **b. Volcano plots:** log2 fold change (reproductive / non-reproductive) against −log10 FDR-adjusted P. Red, padj < 0.05 (dashed line).
- **c. Orthogroup overlaps:** shared caste-biased orthogroups for each pair of species (Excon run, OrthoFinder v3.1.4). Stars: one-sided hypergeometric P (* < 0.05, ** < 0.01, *** < 0.001), universe = orthogroups with a tested gene in both species.
