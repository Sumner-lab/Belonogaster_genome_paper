# Belonogaster_genome_paper

Code used in the *Belonogaster juncea* genome paper. Each folder has its own README with inputs, commands and results.

| Folder | What it covers | Paper |
|---|---|---|
| `Genome assembly and annotation/` | Assembly (Juicer / 3D-DNA), annotation of *B. juncea*, re-annotation of *L. flavolineata* and *P. dominula* | Methods |
| `Mitochondria/` | MitoHiFi run | Methods |
| `HiC_and_satellite/` | Omni-C contact map, Figure 1 assembly, the 296-bp satellite and the search for it in 11 other vespids | Figure 1; Results, *Transposable elements and repeats* |
| `TE_checks/` | RepeatMasker checks: how much of the satellite RepeatMasker masks, and transposable-element overlap of the *B. juncea*-amplified gene families | Results, *Transposable elements and repeats* and *Gene expansions and contractions* |
| `CAFE_Figure4/` | Plotting code and tables for the gene family expansion / contraction figure (Excon v2.5.1, CAFE5) | Figure 4, Supplementary Figure 2 |
| `bin/` | Chromosome GO enrichment summary (`summarize_chromosome_go.3.R`) and `ChopGO.pl` | Figure 5 |
| `GeneOntologyOverlap/` | GO enrichment of caste-biased genes and its overlap between species | Results, RNA-seq |
| `DEG_overlap/` | Overlap of caste-biased orthogroups between species (hypergeometric tests) | Figure 6c |
| `Figure6/` | DESeq2 refit from the raw counts, PCA and volcano plots, assembly of Figure 6 | Figure 6 |

Large inputs (assembly, annotations, CAFE and GO results, genome reports) are on NCBI and Zenodo; the DOIs are given in each README and in the paper's Data availability section.
