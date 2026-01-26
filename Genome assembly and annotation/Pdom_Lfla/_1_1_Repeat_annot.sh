#!/bin/bash

### Setting variables and directories ------------------------------------------------------------

SPECIES=$1
GENOME=$2
GFF=$3
THREADS=$4
JOB_NAME=$5

STEP=02_TE_annotation
prefix=""

BASENAME_GENOME=$(basename ${GENOME/.fna.gz/})

WORKDIR=$HOME/DATA
BIN=$HOME/scripts/bin/
TE_DIR=$WORKDIR/$STEP$prefix/1_HiTE
TEMP=$WORKDIR/tmp/$STEP
SCRIPTS=$HOME/scripts/02_TE_annotation/

for dirs in $TE_DIR/ $TEMP/
do
  if [ ! -d "$dirs" ]
  then
	mkdir -p $dirs
  fi
done

#### Options -------------------------------------------------------------------------------------

DB_PATH=$WORKDIR/databases/
TE_LIB=$DB_PATH/dfam/dfam38_repeatmasker.metazoa.curated.fa

PASSED=$WORKDIR/01_QC$prefix/busco_passed/busco_genome_passed.csv

HITE_PARAM=$(echo "--thread $THREADS --plant 0 --recover 1 --annotate 1")

#### Modules and envs ----------------------------------------------------------------------------
MODULES_LOAD(){
to_load=$1

. $BIN/MODULES_LOAD.sh $to_load
}
#### Starting ------------------------------------------------------------------------------------
# CLEANUP
TPM_CLEAN(){
rm -r $TEMP/
}
# Unzipping
UNZ_GFFFASTA(){
GFF_FASTA=$1

. $BIN/UNZ_GFFFASTA.sh ${GENOME} ${GFF} ${GFF_FASTA}
}
# HiTe on genomes
TE_ANNOT(){
MODULES_LOAD HITE
HITE_MAIN=$CONDA_PREFIX/share/HiTE/main.py

NAME="HiTE"

echo $NAME

UNZ_GFFFASTA

if [ -z "$(grep $(basename $GENOME) $PASSED)" ]
then
echo "WARNING: $(basename $INPUT_GEN) didn't pass thresholds."
exit
fi

if [ ! -f "$TE_DIR/$(basename ${INPUT_GEN/.fna/})/HiTE.tbl" ] || [ ! -s "$TE_DIR/$(basename ${INPUT_GEN/.fna/})/HiTE.tbl" ] || [ ! -f "$TE_DIR/$(basename ${INPUT_GEN/.fna/})/HiTE.gff" ] || [ ! -s "$TE_DIR/$(basename ${INPUT_GEN/.fna/})/HiTE.gff" ]
then

  echo $'\n'"Running HiTe on genomes at $(date)."

  UNZ_GFFFASTA FASTA

  python \
  $HITE_MAIN \
  --genome $INPUT_GEN \
  --outdir $TE_DIR/$(basename ${INPUT_GEN/.fna/}) \
  --curated_lib $TE_LIB \
  $HITE_PARAM

  echo $'\n'"HiTE on genomes finished at $(date)."
else
  echo "Results already exist."
fi
}
#### COMMANDS
TE_ANNOT
TPM_CLEAN
