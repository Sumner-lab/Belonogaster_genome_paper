#!/bin/bash

GENOME=$1
GFF=$2
GFF_FASTA=$3

export INPUT_GEN=$TEMP/$(basename ${GENOME/.fna.gz/.fna})
export INPUT_GFF=$TEMP/$(basename ${GFF/.gff.gz/.gff})

if [ "$GFF_FASTA" == "BOTH" ]
then
if [ ! -f "$INPUT_GEN" ] || [ ! -s "$INPUT_GEN" ]
then
gunzip -c $GENOME > $INPUT_GEN
fi
if [ ! -f "$INPUT_GFF" ] || [ ! -s "$INPUT_GFF" ]
then
gunzip -c $GFF > $INPUT_GFF
fi
elif [ "$GFF_FASTA" == "GFF" ]
then
if [ ! -f "$INPUT_GFF" ] || [ ! -s "$INPUT_GFF" ]
then
gunzip -c $GFF > $INPUT_GFF
fi
elif [ "$GFF_FASTA" == "FASTA" ]
then
if [ ! -f "$INPUT_GEN" ] || [ ! -s "$INPUT_GEN" ]
then
gunzip -c $GENOME > $INPUT_GEN
fi
fi

