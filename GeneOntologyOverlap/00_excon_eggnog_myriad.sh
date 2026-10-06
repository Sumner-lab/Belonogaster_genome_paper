#!/bin/bash
## ============================================================
## GO annotation of the final gene models with eggNOG-mapper, via excon
## Run on UCL Myriad, 2026-10-05, in GOTITS_excon/excon_caste_go/excon
## (a copy of Eco-Flow/excon at v2.5.0-15-g3d77413; commits after the v2.5.0 tag
## change only CI, tests and README).
##
## Outputs used downstream (copied to go_files_excon_v2.5.0/):
##   <species>.go.txt                gene-level GO (EGGNOG_TO_GO)
##   <species>.emapper.annotations   eggNOG-mapper 2.1.13 output
## ============================================================

## input.csv: species, genome FASTA, final annotation GFF (no header).
## These are the same files as in Input_genomes_and_gff/ (GFF md5, uncompressed):
##   Belonogaster_juncea_final_annot.gff          cd864b470d5be5feabbd0d83b8f3755f
##   Liostenogaster_flavolineata_final_annot.gff  2d943adb036cd8e55fe496977801640c
##   Polistes_dominula_final_annot.gff            c25f14a90d5ea81c6b7dd8574abce42f
## On Myriad the lines were taken from an earlier run that used these GFFs:
##   grep -h -E "^(Belonogaster_juncea|Liostenogaster_flavolineata|Polistes_dominula)," \
##     ../excon_belonogaster3/*.csv | sort -u > input.csv
## Template (replace paths):
##   Belonogaster_juncea,/path/B_junceaGenome.FINAL.fasta.gz,/path/Belonogaster_juncea_final_annot.gff
##   Liostenogaster_flavolineata,/path/SoftmaskedL_flavolineata_Dovetail.FINAL.fasta.gz,/path/Liostenogaster_flavolineata_final_annot.gff
##   Polistes_dominula,/path/P_dominulaGenome.FINAL.fasta.gz,/path/Polistes_dominula_final_annot.gff
##
## Note: an earlier excon run used L_flavolineataBraker.gff.gz (16,112 genes) and
## P_dominula.gff.gz (12,525 genes). Those are different gene models from the
## ones used for the counts, and gave GO files whose gene IDs did not match.
## Always check new GO files with 00_validate_go_annotation.py.

## myriad.config is the local cluster config; --custom_config is resolved
## relative to excon's nextflow.config, so keep it inside the excon folder.
nextflow run main.nf -profile singularity -bg \
  --custom_config myriad.config \
  --input input.csv --outdir results_caste_go \
  --run_eggnog --eggnog_data_dir /home/ucbtcdr/Scratch/GOTITS_excon/excon_updates_17Mar26/eggnog_data \
  --eggnog_tax_scope 50557 --eggnog_score 60 --eggnog_pident 40 \
  --eggnog_query_cover 40 --eggnog_subject_cover 40 \
  --skip_cafe

## Collect the GO files
mkdir -p go_for_mac
find results_caste_go \( -name "*.go.txt" -o -name "*.emapper.annotations" \) -exec cp {} go_for_mac/ \;
