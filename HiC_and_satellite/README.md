# Omni-C contact map (Figure 1b) and the 296-bp satellite

This folder holds the code behind three results:
- the contact map of the final *B. juncea* assembly (Figure 1b);
- the satellite DNA results (Tandem Repeats Finder);
- the search for the satellite in other vespid genomes.

## Files

| File | What it does |
|---|---|
| `replot_omnic_final.py` | Figure 1b: Omni-C contact map of chromosomes 1–39 in the final assembly |
| `run_trf.sh` | Runs Tandem Repeats Finder on the final assembly (UCL Myriad job script) |
| `trf_summary.py` | Summarises the TRF output: satellite content, the 296-bp family, overlap with Omni-C contacts and genes |
| `satellite_scan.py` | Looks for the satellite in other genomes |
| `satellite_296bp_consensus.fa` | Consensus monomer of the satellite (296 bp), written by `trf_summary.py` |

## Inputs

Put these files in this folder before running the scripts.

| File | Source |
|---|---|
| `B_junceaGenome.FINAL.fasta` | Final assembly, NCBI GCA_057814355.1 |
| `Belonogaster_juncea_final_annot.gff.gz` | Annotation V2, Zenodo 10.5281/zenodo.23169178 |
| `B_junceaGenome.0.hic`, `B_junceaGenome.0.assembly`, `B_junceaGenome.0.review.assembly` | Juicer/3D-DNA contact map and layouts from the Juicebox curation |
| `N0.tsv` | Orthogroups from the excon run, in the CAFE files on Zenodo; used only for the gene overlap in `trf_summary.py` |
| `B_juncea_trf.txt` | Output of `run_trf.sh` |

## 1. Contact map (Figure 1b)

```
pip install hic-straw numpy matplotlib
python replot_omnic_final.py
```

The `.hic` file is drawn in the original 3D-DNA layout (`.0.assembly`). The script moves each 50-kb bin to its place in the curated layout (`.0.review.assembly`). It handles fragments that were split or inverted during curation.

It then matches the curated scaffolds to the final `HiC_scaffold_N` names by length. 3D-DNA adds a 500-bp gap at each join, which the matching allows for.

Finally it sums the contacts into 250-kb bins and plots chromosomes 1–39 in numerical order. It writes three versions:
- `Figure1b_OmniC_final_assembly`: log scale, all contacts;
- `…_observed`: observed counts; bin pairs below the 75th percentile of between-chromosome contacts are white;
- `…_clean`: balanced by iterative correction; contacts up to the 99th percentile of between-chromosome contacts are white.

It also prints the share of contacts within chromosomes and the contact density per chromosome.

## 2. Tandem repeats and the 296-bp satellite

```
qsub run_trf.sh                    # TRF v4.09.1, ~a few hours with 8 cores
python trf_summary.py              # needs replot_omnic_final.py in this folder (reuses its layout code)
```

`trf_summary.py` merges overlapping arrays and reports the tandem-repeat content of chromosomes 1–39 by class:
- microsatellite (< 10 bp);
- minisatellite (10–99 bp);
- the 296-bp satellite: repeat units of 285–305 bp, or two or three times that length;
- other satellites.

It also checks three more things.
- **One family?** It counts how many long arrays (≥ 10 kb) share at least half of their 10-mers with the consensus monomer of the longest array. Both strands and all rotations are compared.
- **Hi-C overlap:** where the satellite sits relative to the Omni-C contacts, in 50-kb windows.
- **Genes:** how many genes lie inside, or within 10 kb of, satellite arrays.

Results for the final assembly:
- the satellite covers 63.7 Mb, 20.5% of chromosomes 1–39, in 685 arrays on all 39 chromosomes; the longest is 2.79 Mb;
- 91% of long arrays belong to one family;
- 87% of the satellite lies in windows with < 10% of the median number of Omni-C contacts.

## 3. Is the satellite found in other genomes?

```
python satellite_scan.py satellite_296bp_consensus.fa genome1.fa.gz genome2.fa.gz ...
```

The script takes every 12-mer of the monomer, on both strands and from the circular sequence. For each genome it reports how much sequence falls in 2-kb windows where more than 1% (or 10%) of positions match one of these 12-mers.

Random sequence gives about 0.004%. A homologous satellite array gives several percent, even at roughly 75–80% identity.

| Genome | Windows > 1% | Densest window |
|---|---|---|
| *B. juncea* | 65.9 Mb | 85.7% |
| *P. dominula* | 0 | 0.7% |
| *L. flavolineata* | 0 | 0.9% |

## Software

- Python 3.11 with numpy, pandas, matplotlib, and hic-straw 1.3.1 (for `.hic` files).
- Tandem Repeats Finder v4.09.1 (Benson 1999).
