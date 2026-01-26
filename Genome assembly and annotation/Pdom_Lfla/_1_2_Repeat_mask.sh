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
TE_DIR=$WORKDIR/$STEP$prefix
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

SMALL_TE_SIZE=1000 # max size of small TE which will not be masked

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
NAME="MASK"

echo $NAME
UNZ_GFFFASTA

for DIR_OUT in 1_HiTE #2_RM
do

MODULES_LOAD BEDTOOLS

if [ ! -d "$TE_DIR/$DIR_OUT/masked" ]
then
	mkdir -p $TE_DIR/$DIR_OUT/masked
fi

if [ "$DIR_OUT" == "1_HiTE" ]
then
OUT_TBL="HiTE.tbl"
OUT_GFF="HiTE.gff"
elif [ "$DIR_OUT" == "2_RM" ]
then
OUT_TBL="RM.tbl"
OUT_GFF="RM.gff"
fi


SUM_GEN(){
echo "Length hardmasked: $(grep -v ">" $TE_DIR/$DIR_OUT/masked/$(basename $INPUT_GEN) | grep -o "N" | wc -c)" \
> $TE_DIR/$DIR_OUT/masked/$(basename ${INPUT_GEN/.fna/}.masking_stats.txt)

echo "Length softmasked: $(grep -v ">" $TE_DIR/$DIR_OUT/masked/$(basename $INPUT_GEN) | grep -o "a\|c\|t\|g" | wc -c)" \
>> $TE_DIR/$DIR_OUT/masked/$(basename ${INPUT_GEN/.fna/}.masking_stats.txt)

echo "Length not masked: $(grep -v ">" $TE_DIR/$DIR_OUT/masked/$(basename $INPUT_GEN) | wc -c)" \
>> $TE_DIR/$DIR_OUT/masked/$(basename ${INPUT_GEN/.fna/}.masking_stats.txt)

echo "Percentage of the genome masked: $(echo "$(grep -v ">" $TE_DIR/$DIR_OUT/masked/$(basename $INPUT_GEN) | grep -o "N" | wc -c) * 100 / $(grep -v ">" $TE_DIR/masked/$(basename $INPUT_GEN) | wc -c)" | R --no-save | grep "\[1\]")" \
>> $TE_DIR/$DIR_OUT/masked/$(basename ${INPUT_GEN/.fna/}.masking_stats.txt)
}

if [ -f "$TE_DIR/$DIR_OUT/$(basename ${INPUT_GEN/.fna/})/$OUT_GFF" ]
then
  echo $'\n'"Masking genome at $(date)."

  if [ "$DIR_OUT" == "1_HiTE" ]
  then
    UNZ_GFFFASTA FASTA

    grep -v "A-rich\|C-rich\|G-rich\|T-rich\|tRNA\|rRNA\|)n" \
    $TE_DIR/$DIR_OUT/$(basename ${INPUT_GEN/.fna/})/$OUT_GFF | \
	awk -v TE_small="$SMALL_TE_SIZE" 'BEGIN{FS=OFS="\t"} {if ( $5 - $4 > TE_small) {print $0}}' \
    > $TE_DIR/$DIR_OUT/$(basename ${INPUT_GEN/.fna/})/${OUT_GFF}_filtr

    if [ ! -f "$TE_DIR/$DIR_OUT/masked/$(basename $INPUT_GEN)" ] || [ ! -s "$TE_DIR/$DIR_OUT/masked/$(basename $INPUT_GEN)" ]
    then

      bedtools maskfasta -soft -fi $INPUT_GEN \
      -bed $TE_DIR/$DIR_OUT/$(basename ${INPUT_GEN/.fna/})/${OUT_GFF}_filtr \
      -fo $TE_DIR/$DIR_OUT/masked/$(basename $INPUT_GEN)


      MODULES_LOAD R
      SUM_GEN

    fi
  elif [ "$DIR_OUT" == "2_RM" ]
  then

    cp $TE_DIR/$DIR_OUT/$(basename ${INPUT_GEN/.fna/})/$(basename ${INPUT_GEN}).masked \
	$TE_DIR/$DIR_OUT/masked/$(basename $INPUT_GEN)

	MODULES_LOAD R
    SUM_GEN
  fi
  echo $'\n'"Masking finished at $(date)."
else
  echo "$OUT_GFF doesn't exist."
fi

done
}
#### COMMANDS
TE_ANNOT
TPM_CLEAN
