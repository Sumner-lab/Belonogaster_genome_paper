# Figure 4 and Supplementary Figure 2: gene family expansions and contractions

The CAFE5 analysis was run with Excon v2.5.1 (https://github.com/Eco-Flow/excon, commit b809d7c) on the 12 species, using OrthoFinder v3.1.4 orthogroups. The selected model was four gamma rate categories with a Poisson root distribution. Full results are on Zenodo: https://doi.org/10.5281/zenodo.23373652.

This folder holds the plotting scripts and the small tables they read.

## Files

| File | What it is |
|---|---|
| `plot_go_tree_heatmap.R` | Species tree with the number of significant expansions (+) and contractions (−) on each branch, next to a heatmap of enriched GO terms (Figure 4a, b; Supplementary Figure 2) |
| `plot_family_tree.R` | Copy number of one gene family at each tip and ancestral node, with significant gains and losses (Figure 4c) |
| `SpeciesTree_rooted_node_labels.txt` | OrthoFinder species tree (`orthofinder_cafe/ortho_cafe/Species_Tree/`) |
| `cafe_numbered_tree.nwk` | The same tree with internal nodes labelled with their CAFE ids (from CAFE's `Gamma_asr.tre`) |
| `CAFE_summary.txt` | Significant expansions / contractions per branch (Supplementary Table 1: CAFE) |
| `GO_expansion_terms_BH.tsv`, `GO_contraction_terms_BH.tsv` | Benjamini–Hochberg-corrected P of each GO term (rows) in each species (columns), from the topGO tables in the Excon `cafe_go` output; NA = not nominally significant. Used in the paper |
| `GO_*_terms_bonferroni.tsv` | The same with Bonferroni correction |

## Commands

R 4.4 with ape, ggtree, ggplot2 and patchwork. A UTF-8 locale is needed for the − sign.

Figure 4a, b:
```
LC_ALL=en_GB.UTF-8 Rscript plot_go_tree_heatmap.R SpeciesTree_rooted_node_labels.txt GO_expansion_terms_BH.tsv \
  Figure4ab_CAFE_expansions_BH CAFE_summary.txt --cafe-tree cafe_numbered_tree.nwk --pval 0.05 --min-species 1 --top all \
  --fill-label "-log10(P)\nBH (FDR)" --tip-size 6 --label-size 13 --ann-size 5 --tree-frac 2 --heat-frac 2.2 --label-pad 0.75 \
  --cladogram --grey-branches "Liostenogaster_flavolineata,14" --caption none
```
Supplementary Figure 2: the same command with `GO_contraction_terms_BH.tsv` and the prefix `SuppFigure2_CAFE_contractions_BH`.

Figure 4c (largest olfactory receptor family, OG0000009). `Out_cafe_k4_poisson` is the CAFE5 output folder of the selected model:
```
LC_ALL=en_GB.UTF-8 Rscript plot_family_tree.R SpeciesTree_rooted_node_labels.txt cafe_numbered_tree.nwk \
  Out_cafe_k4_poisson OG0000009 Figure4c_OG0000009_final "Liostenogaster_flavolineata,14,1"
```
The panels were joined with ImageMagick.

## Caveat: root branches

*L. flavolineata* is the only species on one side of the root, so changes on its branch and on the branch leading to the other 11 species (node 14) cannot be separated, and their counts depend on the model:

| Model | *L. flavolineata* | Node 14 |
|---|---|---|
| Selected model | 1,436 up / 0 down | 117 / 0 |
| Base model | 1,373 / 201 | 10 / 161 |

These two branches are greyed in the figure and not interpreted.
