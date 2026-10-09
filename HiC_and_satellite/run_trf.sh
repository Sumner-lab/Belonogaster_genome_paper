#!/bin/bash -l
# Tandem Repeats Finder on the final B. juncea assembly (UCL Myriad, SGE). Submit with: qsub run_trf.sh
# TRF v4.09.1 binary: https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
#   (saved as ./trf and made executable). Parameters: match 2, mismatch 7, indel 7, PM 80, PI 10,
#   min score 50, max period 2000 bp; -l 10 allows arrays up to 10 Mb; -ngs writes one compact table.
#$ -N trf_bjuncea
#$ -l h_rt=24:00:00
#$ -l mem=4G
#$ -pe smp 8
#$ -cwd

GENOME=B_junceaGenome.FINAL.fasta.gz        # final assembly (NCBI GCA_057814355.1)

zcat "$GENOME" > genome.fa
awk '/^>/{f=substr($1,2)".fa"} {print > f}' genome.fa          # one file per scaffold, run 8 at a time
ls HiC_scaffold_*.fa | xargs -P 8 -I{} sh -c './trf {} 2 7 7 80 10 50 2000 -h -ngs -l 10 > {}.trf.txt'
cat HiC_scaffold_*.fa.trf.txt > B_juncea_trf.txt
