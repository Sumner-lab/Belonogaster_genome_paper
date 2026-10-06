# GeneOntologyOverlap

GO enrichment of caste-biased genes in *Belonogaster juncea*, *Liostenogaster flavolineata* and *Polistes dominula*, and the overlap of enriched GO terms between species. Produces Figure 7, Supplementary Table 3 and Supplementary Figure 3.

## Scripts

| Script | Runs on | What it does |
|---|---|---|
| `00_excon_eggnog_myriad.sh` | HPC (UCL Myriad) | Record of how the GO files were made: [excon](https://github.com/Eco-Flow/excon) v2.5.0 with eggNOG-mapper v2.1.13 on the final gene annotations |
| `00_validate_go_annotation.py` | local | Checks that the gene IDs in the GO files match the gene IDs in the DESeq2 results |
| `01_run_topGO.R` | local | Builds the caste-biased gene lists, runs GO enrichment, writes Supplementary Table 3 |
| `02_plot_GO_figure.R` | local | Figure 7 (dot plot of top enriched terms) |
| `03_go_overlap.py` | local | Overlap of enriched terms between gene sets, permutation test, Supplementary Figure 3 |

## Inputs

Paths are set at the top of each script (`CONFIG`/`ORTHOGROUPS` in `00_validate_go_annotation.py`, `SPECIES` in `01_run_topGO.R`). Edit them to point at your copies.

- **GO annotations**, one gene-level file per species (two columns, gene ID and GO ID, tab-separated), from Zenodo: [10.5281/zenodo.20608364](https://doi.org/10.5281/zenodo.20608364). For *L. flavolineata* and *P. dominula*, use the files made from the final gene annotations (added in version 2 of the record).
- **DESeq2 results**, one table per species with columns `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj` (Supplementary Table 2). Positive log2 fold change = higher in reproductives.
- **Final gene annotations** (`*_final_annot.gff`), used only by `00_validate_go_annotation.py`.
- **OrthoFinder `Orthogroups.tsv`**, used only by `00_validate_go_annotation.py`.

## Requirements

- R ≥ 4.4 with topGO, writexl, ggplot2, dplyr, forcats, stringr and ggh4x
- Python 3:
  - `00_validate_go_annotation.py` needs only the standard library
  - `03_go_overlap.py` needs pandas, openpyxl, matplotlib, numpy and scipy

## Run

```
python3 00_validate_go_annotation.py
Rscript 01_run_topGO.R
Rscript 02_plot_GO_figure.R
python3 03_go_overlap.py Supplementary_Table_3_GO_enrichment.xlsx -o go_overlap --alpha 0.01 \
    --manifest results/overlap_manifest.tsv
```

`03_go_overlap.py` runs 10,000 permutations by default (seed 1, about 2 minutes). Without `--manifest`, it only describes the overlap (summary tables and UpSet plot) and runs no tests.

## Method summary

**Gene sets**
- DEGs have DESeq2 adjusted P < 0.05.
- Queen-biased = log2 fold change > 0; worker-biased = log2 fold change < 0.

**Enrichment**
- topGO 'classic' algorithm with one-sided Fisher's exact tests, run separately for each species, caste and GO domain (BP, MF, CC).
- Background: expressed genes, i.e. genes in the DESeq2 results table, that have at least one GO annotation.
- The same analysis with all annotated genes as the background is reported as a robustness check (`results/summary.tsv`).
- P-values are Benjamini–Hochberg corrected across all terms tested within each domain.

**Figure 7**
- Shows the 12 most significant terms (P < 0.01) per species and caste.
- Black outline = FDR < 0.05.

**Overlap test**
- Uses terms with P < 0.01, and only pairs of gene sets from different species.
- GO terms are not independent, so the null distribution comes from permutation. In each permutation:
  - random gene sets of the observed sizes are drawn from each species' background (queen- and worker-biased sets without overlap);
  - all terms are re-tested;
  - the same number of top-ranked terms as observed is kept.
- P-values are Benjamini–Hochberg corrected across the between-species pairs.
- A hypergeometric test is also reported for reference. It treats terms as independent, so it is anti-conservative.

## Outputs

- `gene_lists/`: queen- and worker-biased genes per species
- `results/`: all tested terms per gene set (`*_all_terms.tsv`), `summary.tsv`, `annotation_validation.tsv`, gene-to-term membership and `overlap_manifest.tsv` (inputs for step 3)
- `Supplementary_Table_3_GO_enrichment.xlsx`: terms with P < 0.05 per gene set, with the DEGs in each term
- `Figure7_GO_dotplot_caste.pdf`
- `go_overlap/`:
  - `pairwise_overlap.tsv`: test results per pair
  - `overlap_figure.pdf`: Supplementary Figure 3
  - `upset.pdf` and `overlap_stats.pdf`: the two panels as separate files
  - shared-term tables and the permutation null distributions
