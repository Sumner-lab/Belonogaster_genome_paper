# RepeatMasker checks on the *B. juncea* genome

The RepeatMasker run is part of nf-core/genomeqc v1.0.1 (`--te repeatmasker --repeatmasker_speed qq --famdb_lineage hymenoptera`, Dfam 3.9 partitions 0 and 3). It uses RepeatMasker 4.1.5 with the Dfam families for Hymenoptera (curated and uncurated) and no de novo library.

## Inputs (put them in this folder)

| File | Source |
|---|---|
| `Belonogaster_juncea.out` | RepeatMasker output from genomeqc (`repeatmasker/Belonogaster_juncea/`); genome reports on Zenodo, https://doi.org/10.5281/zenodo.23378133 |
| `B_juncea_trf.txt` | Tandem Repeats Finder output, from `HiC_and_satellite/run_trf.sh` |
| `Belonogaster_juncea_final_annot.gff` | Annotation V2, https://doi.org/10.5281/zenodo.23169178 |
| `N0_OFv2.tsv` | Hierarchical orthogroups (N0) from OrthoFinder v2.5.5, from the Excon robustness run |

## Scripts

```
python rm_check.py     # how much of the 296-bp satellite RepeatMasker masks
python amp_te.py       # TE overlap of the coding sequence of amplified gene families
```

## Results

- RepeatMasker masks 32.0% of the *B. juncea* genome: DNA 7.2%, LTR 3.4%, LINE 2.1%, unclassified 7.3%, RC 0.4%, simple repeats 10.1%, low complexity 1.3%, satellites 0.09%.
- The 296-bp satellite covers 63.8 Mb of chromosomes 1–39, but **RepeatMasker masks only 2.7% of it**, so it is not in Dfam. RepeatMasker plus the satellite comes to about 52% of the genome, close to the de novo estimate from HiTE (50.3%).
- Amplified families: N0 orthogroups where *B. juncea* has ≥ 20 copies and ≥ 5 times as many as any other species (34 orthogroups, 1,797 genes).

| Gene set | Median CDS covered by TEs | Genes ≥ 50% covered | Genes with any overlap |
|---|---|---|---|
| Amplified | 0.89 | 57% | 65% |
| Single-copy in all 12 species (4,593) | 0.00 | 3% | 40% |
| Random 3,000 genes | 0.00 | 19% | 48% |

The TE classes in the amplified genes are mainly LINE, unclassified, DNA and LTR.

## Caveat: Dfam coverage

Dfam 3.9 has species-specific families for 10 of the 12 species (*A. nigricornis* 1,378; *P. dominula* 745; *V. germanica* 692; *V. vulgaris* 628; *D. sylvestris* 590; *D. saxonica* 553; *D. media* 504; *V. pensylvanica* 295; *V. crabro* 277; *V. velutina* 260), but none for *B. juncea* or *L. flavolineata*. RepeatMasker therefore probably under-estimates repeats in these two species relative to the others.
