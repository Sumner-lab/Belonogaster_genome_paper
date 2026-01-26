#!/bin/bash

### Setting variables and directories ------------------------------------------------------------

SPECIES=$1
GENOME=$2
GFF=$3
THREADS=$4
JOB_NAME=$5

STEP=03_Annotation
prefix=""

BASENAME_GENOME=$(basename ${GENOME/.fna.gz/})

CLADE=hymenoptera
TE_SOFT="HiTE" #HiTE or RM

WORKDIR=$HOME/DATA
BIN=$HOME/scripts/bin/
QC_DIR=$WORKDIR/01_QC$prefix
ANNOT_DIR=$WORKDIR/$STEP$prefix/RNAprot_${CLADE}/ 
TE_AN_DIR=$WORKDIR/02_TE_annotation$prefix
TEMP=$WORKDIR/tmp
TRANCRIPTOMEFILE_NCBI=$WORKDIR/$STEP$prefix/transcriptomefile_forannot.tsv
TRANCRIPTOMEFILE_RDS=$WORKDIR/00_DATA$prefix/datasets_RDS_download/transcriptome_list.txt
SCRIPTS=$HOME/scripts/03_Genome_annotation


for DIRS in $ANNOT_DIR/transcriptomes/raw_data $ANNOT_DIR/long_iso/gff $ANNOT_DIR/long_iso/interpro $ANNOT_DIR/long_iso/passed $ANNOT_DIR/busco_proteome $ANNOT_DIR/gff_stats
do
if [ ! -d "$DIRS" ]
then
        mkdir -p $DIRS
fi
done

if [ $TE_SOFT == "HiTE" ]
then
  TE_DIR=$TE_AN_DIR/1_HiTE/
elif [ $TE_SOFT == "RM" ]
then
  TE_DIR=$TE_AN_DIR/2_RM/
fi

#### Options -------------------------------------------------------------------------------------

DB=$WORKDIR/databases
DB_NAME=hymenoptera
DB_V=odb10

AUGUSTUS_CONFIG_DIR_INIT=$WORKDIR/envs/conda/braker3_3.0.8/config
AUGUSTUS_SPECIES=$AUGUSTUS_CONFIG_DIR/species/

PROT_DB=$DB/hymenoptera_genomes/Hymenoptera_protein_COMPLETE.faa

RNAMAPPING_CONFIG_PATH=$HOME/config_files/RNA_filtering_mapping.cfg

export BLASTDB=$NT_DB:$NR_DB

READ_THRESHOLD=5000000 #min nb of reads required in bam files

INTERPROSCAN=interproscan.sh
COMPLEASM=$WORKDIR/envs/conda/braker3_3.0.8/bin/compleasm

FILTER_THRESHOLD=( 7 5 5 ) # missing fragmented duplicated

BRAKER_PARAM=$(echo "--threads $THREADS --gff3 --useexisting")
BEST_COMPLEASM_PARAM=$(echo "-n 4")
INTERPROSCAN_PARAM=$(echo "--goterms -ms 300 -T $TEMP")
BUSCO_PARAM_P=$(echo "-m prot -c $THREADS -f -e 1e-05 -l $DB_NAME --datasets_version $DB_V --download_path $DB/busco --offline")

#### Modules and envs ----------------------------------------------------------------------------
MODULES_LOAD(){
to_load=$1
do_stack=$2

. $BIN/MODULES_LOAD.sh $to_load $do_stack
}
#### Starting ------------------------------------------------------------------------------------
# CLEANUP
TPM_CLEAN(){
rm $TEMP/$(basename ${GENOME/.fna.gz/})*
rm $TEMP/scripts/sbatch_${JOB_NAME}*.${BASENAME_GENOME}.*.sh
}
# Unzipping
UNZ_GFFFASTA(){
GFF_FASTA=$1

. $BIN/UNZ_GFFFASTA.sh ${GENOME} ${GFF} ${GFF_FASTA}
}

RNA_DL(){
RESERVE=$1

echo "Downloading RNA seq data for $SPECIES at $(date)"

if [ ! -d "$ANNOT_DIR/transcriptomes/raw_data/${BASENAME_GENOME}/" ]
then
  mkdir -p $ANNOT_DIR/transcriptomes/raw_data/${BASENAME_GENOME}/
fi

#. $SCRIPTS/_download_RDS_data.sh $SPECIES $GENOME $GFF $THREADS $ANNOT_DIR $TRANCRIPTOMEFILE_RDS ## download from UCL internal storage space

. $SCRIPTS/_download_SRA.sh $SPECIES $GENOME $GFF $THREADS $ANNOT_DIR $TRANCRIPTOMEFILE_NCBI $RESERVE

if [ -z "$(ls $ANNOT_DIR/transcriptomes/raw_data/${BASENAME_GENOME}/*.fastq*)" ]
then
echo "something failed: no RNAseq fastq file available..."
exit
fi
}
# Mapping RNA
RNA_MAPPING(){
ROUND=$1

if [ "$ROUND" -eq 1 ]
then
  TAKE_RESERVE="no"
else
  TAKE_RESERVE="yes"
fi

RNA_DL $TAKE_RESERVE

if [ "$(ls $ANNOT_DIR/transcriptomes/raw_data/${BASENAME_GENOME} | wc -l)" == 0 ]
then
  echo "no RNA data available for ${BASENAME_GENOME}"
  exit
fi

PIPE1_OUT=01_trimmed_filtered_data
PIPE2_OUT=02_mapping_assembly

UNZ_GFFFASTA FASTA

sed "s+^WORK_DIR:.*#+WORK_DIR: $ANNOT_DIR/transcriptomes/ #+;s+^RAW_DATA_DIR:.*#+RAW_DATA_DIR: $ANNOT_DIR/transcriptomes/raw_data/ #+;s+^GENOME_FILE:.*#+GENOME_FILE: $INPUT_GEN #+;s+^N_THREADS:.*#+N_THREADS: $THREADS #+;s+^RAW_DATA_DIR:\(.*\)/ #+RAW_DATA_DIR:\1/${BASENAME_GENOME} #+;s+^BASENAME:.*#+BASENAME: ${BASENAME_GENOME} #+" \
${RNAMAPPING_CONFIG_PATH} \
> $ANNOT_DIR/transcriptomes/${BASENAME_GENOME}_$(basename $RNAMAPPING_CONFIG_PATH)

RNA_MAP_PART1(){
MODULES_LOAD TRIMMOMATIC

$BIN/RNA_mapping_annot.sh \
-c $ANNOT_DIR/transcriptomes/${BASENAME_GENOME}_$(basename $RNAMAPPING_CONFIG_PATH) \
-p 1 -s 1

if [ ! -z "$(ls $ANNOT_DIR/transcriptomes/$PIPE1_OUT/01_trimmomatic/paired*${BASENAME_GENOME}*)" ]
then
 rm -r $ANNOT_DIR/transcriptomes/raw_data/${BASENAME_GENOME}/
fi

MODULES_LOAD FASTQC
MODULES_LOAD MULTIQC stack

$BIN/RNA_mapping_annot.sh \
-c $ANNOT_DIR/transcriptomes/${BASENAME_GENOME}_$(basename $RNAMAPPING_CONFIG_PATH) \
-p 1 -s 2

MODULES_LOAD HISAT
MODULES_LOAD SAMTOOLS stack


$BIN/RNA_mapping_annot.sh \
-c $ANNOT_DIR/transcriptomes/${BASENAME_GENOME}_$(basename $RNAMAPPING_CONFIG_PATH) \
-p 2 -s 1 2> $ANNOT_DIR/transcriptomes/HISAT_ERR_${BASENAME_GENOME}
}
RNA_MAP_PART1

MODULES_LOAD SAMTOOLS

for bams in $(ls $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/*${BASENAME_GENOME}*.bam | xargs -n1 basename )
do
  nb_reads=$(samtools coverage $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/$bams | grep -v "#" | cut -f4 | awk '{sum += $1} END {print sum}')
  echo "nb mapping reads in $bams : $nb_reads"
  if [ "$nb_reads" -lt "$READ_THRESHOLD" ]
  then
    mv $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/$bams \
	$ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/flagged_$bams
  fi
done

nb_bams=$(ls $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/*${BASENAME_GENOME}*.bam | grep -v "flagged_"| xargs -n1 basename | wc -l)
nb_reserve=$(awk -F'\t' -v sp_name="$SPECIES" 'split(sp_name,a,"_") {b=a[1] " " a[2]} $1 ~ b && $2 == "RESERVE" {print $0}' $TRANCRIPTOMEFILE_NCBI | wc -l)

if [ "$nb_bams" -lt 3 ] && [ "$nb_reserve" -gt 0 ] && [ "$TAKE_RESERVE" == "no" ]
then
  RNA_DL yes
  RNA_MAP_PART1
fi


for bams in $(ls $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/*${BASENAME_GENOME}*.bam | grep -v "flagged_" | xargs -n1 basename )
do
  nb_reads=$(samtools view -c -q 255 -F 0x2 $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/$bams | cut -f1 | sort | uniq | wc -l)
  if [ "$nb_reads" -lt "$READ_THRESHOLD" ]
  then
    mv $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/$bams \
	$ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/flagged_$bams
  fi
done

nb_bams=$(ls $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/*${BASENAME_GENOME}*.bam | grep -v "flagged_" | xargs -n1 basename | wc -l)
nb_flagged=$(ls $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/*${BASENAME_GENOME}*.bam | grep "flagged_" | xargs -n1 basename | wc -l)

echo $'\n'"nb of bams available: $(expr $nb_bams + $nb_flagged ), nb of flagged bams : $nb_flagged"$'\n'

if [ "$(expr $nb_bams + $nb_flagged)" -eq 0 ]
then
  echo "no RNA data available for ${BASENAME_GENOME}"
  exit
fi


if [ ! -z "$(ls $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/*${BASENAME_GENOME}*.bam)" ]
then
  rm -r $ANNOT_DIR/transcriptomes/$PIPE1_OUT/01_trimmomatic/*${BASENAME_GENOME}*.log
  rm -r $ANNOT_DIR/transcriptomes/$PIPE1_OUT/01_trimmomatic/*${BASENAME_GENOME}*.gz
fi

MODULES_LOAD SAMTOOLS

$BIN/RNA_mapping_annot.sh \
-c $ANNOT_DIR/transcriptomes/${BASENAME_GENOME}_$(basename $RNAMAPPING_CONFIG_PATH) \
-p 2 -s 2

if [ ! -z "$(ls $ANNOT_DIR/transcriptomes/$PIPE2_OUT/03_sorted_bam_files/*${BASENAME_GENOME}*sorted.bam)" ]
then
  rm -r $ANNOT_DIR/transcriptomes/$PIPE2_OUT/01_hisat_mapping/index/*${BASENAME_GENOME}*
  rm -r $ANNOT_DIR/transcriptomes/$PIPE2_OUT/01_hisat_mapping/*${BASENAME_GENOME}*
  rm -r $ANNOT_DIR/transcriptomes/$PIPE2_OUT/02_bam_files/*${BASENAME_GENOME}*
fi

ls $ANNOT_DIR/transcriptomes/$PIPE2_OUT/03_sorted_bam_files/*${BASENAME_GENOME}*sorted.bam | \
xargs -n1 basename | sed "s/\.bam$//g" | \
tr '\n' ',' | sed "s/,$//" \
> $ANNOT_DIR/transcriptomes/list_bam_${BASENAME_GENOME}.txt

}
# Genome Annotation
# BRAKER
GENOME_ANNOT_BRAKER(){

OUT=$ANNOT_DIR/$BASENAME_GENOME/braker3_out

LAUNCH_BRAKER(){
ATTEMPT=$1

if [ ! -f "$OUT/braker.gtf" ] || [ ! -s "$OUT/braker.gtf" ]
then
  if [ -z "$(ls $ANNOT_DIR/transcriptomes/02_mapping_assembly/03_sorted_bam_files/*${BASENAME_GENOME}*sorted.bam)" ] && [ "$ATTEMPT" -eq 1 ]
  then
    RNA_MAPPING $ATTEMPT
	if [ -z "$(ls $ANNOT_DIR/transcriptomes/02_mapping_assembly/03_sorted_bam_files/*${BASENAME_GENOME}*sorted.bam)" ]
    then
	  echo "Something went wrong: couldn't download or map the transcriptomes"
	  exit
	fi
  fi
fi

  ls $ANNOT_DIR/transcriptomes/02_mapping_assembly/03_sorted_bam_files/*${BASENAME_GENOME}*sorted.bam | \
  xargs -n1 basename | sed "s/\.bam$//g" | \
  tr '\n' ',' | sed "s/,$//" \
  > $ANNOT_DIR/transcriptomes/list_bam_${BASENAME_GENOME}.txt

UNZ_GFFFASTA

if [ ! -d "$ANNOT_DIR/augustus_config/" ]
then
  mkdir -p $ANNOT_DIR/augustus_config
fi

if [ ! -d "$ANNOT_DIR/augustus_config/species/" ]
then
  cp -r $AUGUSTUS_CONFIG_DIR_INIT/* $ANNOT_DIR/augustus_config/
fi

if [ ! -d "$OUT" ]
then
  mkdir -p $OUT
fi

cd $OUT

if [ ! -f "$OUT/braker.gtf" ] || [ ! -s "$OUT/braker.gtf" ]
then
  MODULES_LOAD BRAKER #stack
  AUGUSTUS_CONFIG=$(readlink -f $ANNOT_DIR)/augustus_config/
  AUGUSTUS_SPECIES=$(readlink -f $AUGUSTUS_CONFIG)/species/
  export AUGUSTUS_CONFIG_PATH=$(readlink -f $AUGUSTUS_CONFIG)

  echo "Starting genome annotation using BRAKER at $(date): attempt $ATTEMPT."
  rm -r $AUGUSTUS_SPECIES/${BASENAME_GENOME}*

    if [ ! -f "$PROT_DB" ] 
    then
      echo "$PROT_DB does not exist"
      exit
    fi
 
    if [ ! -f "$TE_DIR/masked/$(basename $INPUT_GEN)" ]
    then
      echo "$TE_DIR/masked/$(basename $INPUT_GEN) does not exist"
      exit
    fi

  if [ ! -f "$OUT/braker/Augustus/augustus.hints.gtf" ] || [ ! -f "$OUT/braker/hintsfile.gff" ] || [ ! -f "$OUT/transcripts_merged.gff" ]
  then

	ls $ANNOT_DIR/transcriptomes/02_mapping_assembly/03_sorted_bam_files/*${BASENAME_GENOME}*sorted.bam || \
	(
	  echo "No transcriptome dataset available for ${BASENAME_GENOME}"
	  exit
	)

	braker.pl --species=${BASENAME_GENOME} \
    --genome=$(readlink -f $TE_DIR)/masked/$(basename $INPUT_GEN) \
    --prot_seq=$(readlink -f $PROT_DB) \
    --rnaseq_sets_ids=$(cat $(readlink -f $ANNOT_DIR)/transcriptomes/list_bam_${BASENAME_GENOME}.txt) \
    --rnaseq_sets_dirs=$(readlink -f $ANNOT_DIR)/transcriptomes/02_mapping_assembly/03_sorted_bam_files/ \
    $BRAKER_PARAM
  fi
fi
}
LAUNCH_BRAKER 1

if [ ! -f "$OUT/final_annot.gff" ] || [ ! -s "$OUT/final_annot.gff" ]
then
  if [ ! -f "$OUT/braker/braker.gtf" ]
  then
    echo "not enough evidence, retrying with transcriptomes in reserve"
    LAUNCH_BRAKER 2
	if [ ! -f "$OUT/braker/braker.gtf" ]
    then
      echo "$OUT/braker/braker.gtf does not exist"
	  #rm -r $ANNOT_DIR/$BASENAME_GENOME/
      exit
	fi
  fi

  if [ ! -f "$OUT/braker.gtf" ]
  then
    cp $OUT/braker/braker.gtf \
    $OUT/braker.gtf
  fi

  cp $OUT/braker/Augustus/augustus.hints.gff3 \
  $OUT/augustus.hints.gff3
  cp $OUT/braker/Augustus/augustus.hints.gtf \
  $OUT/augustus.hints.gtf
  cp $OUT/braker/GeneMark-ETP/genemark.gtf \
  $OUT/genemark.gtf
  cp $OUT/braker/GeneMark-ETP/training.gtf \
  $OUT/training.gtf
  cp $OUT/braker/hintsfile.gff \
  $OUT/hintsfile.gff
  cp $OUT/braker/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff \
  $OUT/transcripts_merged.gff
  
  echo "Fixing GFFs at $(date)"
  
  grep ">" $(readlink -f $TE_DIR)/masked/$(basename $INPUT_GEN) | \
  sed "s/ .*//;s/>//" | awk 'BEGIN{FS=OFS="\t"} {a=gensub("\\|","_","g",$1)} {$2 = a}1' \
  > $OUT/fasta_headers.txt
  
   awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$1;next} {for(i in a)if(index($0,i)) $3=a[i]}1' \
  <(grep -v "#" <(cat  $OUT/braker.gtf $OUT/augustus.hints.gtf $OUT/genemark.gtf $OUT/training.gtf) | cut -f1 | sort | uniq) \
  $OUT/fasta_headers.txt \
  > $OUT/fasta_headers.txt.mod
  
  mv $OUT/fasta_headers.txt.mod \
  $OUT/fasta_headers.txt
  
  for files in augustus.hints.gff3 augustus.hints.gtf hintsfile.gff genemark.gtf training.gtf braker.gtf transcripts_merged.gff
  do
    mv $OUT/$files $OUT/${files}.tmp
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$3]=$1;next} $1 in a{$1=a[$1]}1' \
	$OUT/fasta_headers.txt \
	$OUT/${files}.tmp \
	> $OUT/${files}
  done

  echo "refining filrtering thresholds"

  if [ ! -d "$TEMP/${BASENAME_GENOME}/" ]
  then
    mkdir -p $TEMP/${BASENAME_GENOME}/
  fi
  best_by_compleasm.py -m $(readlink -f $TEMP)/${BASENAME_GENOME}/ \
  -d  $(readlink -f $OUT)/braker/ \
  -g $(readlink -f $TE_DIR)/masked/$(basename $INPUT_GEN) \
  -t $THREADS \
  -p $DB_NAME \
  -c $COMPLEASM \
  $BEST_COMPLEASM_PARAM

  if [ -f "$TEMP/${BASENAME_GENOME}/better.gtf" ] || [ -f "$TEMP/${BASENAME_GENOME}/better.gff" ]
  then
    cp $TEMP/${BASENAME_GENOME}/better.g* \
    $OUT/
  else
    cp $OUT/braker.gtf \
    $OUT/better.gtf
  fi

  if [ ! -f "$OUT/better.gff" ]
  then
  MODULES_LOAD AGAT stack

  agat_convert_sp_gxf2gxf.pl \
  -g $OUT/better.gtf \
  -o $OUT/better.gff
  fi

  if [ ! -f "$WORKDIR/01_QC$prefix/minimap_out/${BASENAME_GENOME}/dup_seq.txt" ] || [ ! -s "$WORKDIR/01_QC$prefix/minimap_out/${BASENAME_GENOME}/dup_seq.txt" ]
  then
    cp $OUT/better.gff \
    $OUT/final_annot.gff
  else
	MODULES_LOAD BEDTOOLS

	grep $'\t'"gene"$'\t' \
	$OUT/better.gff \
	> $OUT/better.gene.gff

	bedtools intersect \
	-wa -f 0.9 \
	-a $OUT/better.gene.gff \
	-b <(cut -f 4-6 $WORKDIR/01_QC$prefix/minimap_out/${BASENAME_GENOME}/dup_seq.txt) | \
	cut -f 9 | sed "s/;.*//;s/ID=//" \
	> $OUT/better.gene_duplicate.txt

	MODULES_LOAD AGAT

	agat_sp_filter_feature_from_kill_list.pl \
	--gff $OUT/better.gff \
	--kill_list $OUT/better.gene_duplicate.txt \
	-o $OUT/final_annot.gff
	fi

  echo "adding UTRs"
  ingenannot -v 2 utr_refine \
  $OUT/final_annot.gff \
  $OUT/transcripts_merged.gff \
  $OUT/final_annot_withUTR.gff

  rm -r $TEMP/${BASENAME_GENOME}/

  for dirs in $(ls -d $OUT/*/ | grep -v "interproscan" | xargs -n1 basename)
  do
    tar czvf $OUT/${dirs}.tar.gz \
    $OUT/${dirs}
    rm -r $OUT/${dirs}
  done

  rm $ANNOT_DIR/transcriptomes/02_mapping_assembly/03_sorted_bam_files/*${BASENAME_GENOME}*sorted.bam
  rm $ANNOT_DIR/transcriptomes/${BASENAME_GENOME}_$(basename $RNAMAPPING_CONFIG_PATH)
  rm -r $ANNOT_DIR/transcriptomes/raw_data/*${BASENAME_GENOME}*
fi
echo "Genome annotation finished at $(date)."
}
# functional annotation
FUNCT_ANNOT(){
METHOD="braker3"

OUT=$ANNOT_DIR/$BASENAME_GENOME/${METHOD}_out

MODULES_LOAD AGAT

UNZ_GFFFASTA

if [ -f "$OUT/final_annot.gff" ] && [ -s "$OUT/final_annot.gff" ]
then

  if [ ! -f "$ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa" ] || [ ! -s "$ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa" ]
  then

    echo $'\n'"Extracting longest isoform at $(date)."

    UNZ_GFFFASTA FASTA

    agat_sp_keep_longest_isoform.pl \
    --gff $OUT/final_annot.gff \
    -o $ANNOT_DIR/long_iso/gff/${BASENAME_GENOME}.long_iso.gff

    agat_sp_extract_sequences.pl \
    --gff $ANNOT_DIR/long_iso/gff/${BASENAME_GENOME}.long_iso.gff \
    -f $INPUT_GEN -p \
    -o $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa

    MODULES_LOAD SEQKIT

    seqkit seq -M 100000 $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa \
    > $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa.mod
    mv $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa.mod \
    $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa

  fi

  if [ ! -f "$ANNOT_DIR/long_iso/interpro/${BASENAME_GENOME}.faa" ]
  then
    sed "s/\*//g" $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa \
    > $ANNOT_DIR/long_iso/interpro/${BASENAME_GENOME}.faa
  fi

  echo $'\n'"Functional annotation at $(date)."

  MODULES_LOAD INTERPROSCAN

  if [ ! -d "$OUT/interproscan/1_interproscan/" ]
  then
    mkdir -p $OUT/interproscan/1_interproscan/
    mkdir -p $OUT/interproscan/2_interproscan/
  fi

  if [ ! -f "$(readlink -f $OUT)/interproscan/2_interproscan/${BASENAME_GENOME}.faa.gff3" ]
  then
    $INTERPROSCAN \
    -appl Pfam -t p -f GFF3 -f TSV \
    -i $(readlink -f $ANNOT_DIR/long_iso/interpro)/${BASENAME_GENOME}.faa \
    -b $(readlink -f $OUT)/interproscan/1_interproscan/ \
    $INTERPROSCAN_PARAM

    $INTERPROSCAN \
    -t p -f GFF3 -f TSV \
    -i $(readlink -f $ANNOT_DIR/long_iso/interpro)/${BASENAME_GENOME}.faa \
    -b $(readlink -f $OUT)/interproscan/2_interproscan/ \
    $INTERPROSCAN_PARAM
  fi
  
  echo "Genome annotation finished at $(date)."

else
  echo "ERROR: $OUT/final_annot.gff doesn't exist"
fi

}
# BUSCO on new proteins
BUSCO_P(){
METHOD="braker3"

OUT=$ANNOT_DIR/$BASENAME_GENOME/${METHOD}_out

MODULES_LOAD AGAT

UNZ_GFFFASTA

if [ -f "$OUT/final_annot.gff" ] && [ -s "$OUT/final_annot.gff" ]
then

  if [ ! -f "$ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa" ] || [ ! -s "$ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa" ]
  then

    echo $'\n'"Extracting longest isoform at $(date)."
    UNZ_GFFFASTA FASTA

    agat_sp_keep_longest_isoform.pl \
    --gff $OUT/final_annot.gff \
    -o $ANNOT_DIR/long_iso/gff/${BASENAME_GENOME}.long_iso.gff

    agat_sp_extract_sequences.pl \
    --gff $ANNOT_DIR/long_iso/gff/${BASENAME_GENOME}.long_iso.gff \
    -f $INPUT_GEN -p \
    -o $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa

    MODULES_LOAD SEQKIT

    seqkit seq -M 100000 $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa \
    > $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa.mod
    mv $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa.mod \
    $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa

  fi

  MODULES_LOAD BUSCO


  if [ ! -f "$ANNOT_DIR/busco_proteome/${BASENAME_GENOME}/run_${DB_NAME}_${DB_V}/short_summary.txt" ]
  then
    echo $'\n'"Running BUSCO on proteomes at $(date)."

    busco -i $ANNOT_DIR/long_iso/${BASENAME_GENOME}.faa \
    -o ${BASENAME_GENOME} \
    --out_path $ANNOT_DIR/busco_proteome/ \
    $BUSCO_PARAM_P

    for files in $(ls -d $ANNOT_DIR/busco_proteome/${BASENAME_GENOME}/run_${DB_NAME}_${DB_V}/*/ | xargs -n1 basename)
    do
      tar czvf \
      $ANNOT_DIR/busco_proteome/${BASENAME_GENOME}/run_${DB_NAME}_${DB_V}/${files}.tar.gz \
      $ANNOT_DIR/busco_proteome/${BASENAME_GENOME}/run_${DB_NAME}_${DB_V}/${files}/
      rm -r $ANNOT_DIR/busco_proteome/${BASENAME_GENOME}/run_${DB_NAME}_${DB_V}/${files}/
    done
  fi


  MISSING_BUSCO=$(grep "C:" $ANNOT_DIR/busco_proteome/${BASENAME_GENOME}/run_${DB_NAME}_${DB_V}/short_summary.txt | sed "s/.*M://;s/%.*//")
  DUPLICATED_BUSCO=$(grep "C:" $ANNOT_DIR/busco_proteome/${BASENAME_GENOME}/run_${DB_NAME}_${DB_V}/short_summary.txt | sed "s/.*D://;s/%.*//")
  FRAGMENTED_BUSCO=$(grep "C:" $ANNOT_DIR/busco_proteome/${BASENAME_GENOME}/run_${DB_NAME}_${DB_V}/short_summary.txt | sed "s/.*F://;s/%.*//")

  if (( $(bc <<<"$MISSING_BUSCO < ${FILTER_THRESHOLD[0]} && $FRAGMENTED_BUSCO < ${FILTER_THRESHOLD[1]} && $DUPLICATED_BUSCO < ${FILTER_THRESHOLD[2]}") ))
  then
    echo ${BASENAME_GENOME} \
    >> $ANNOT_DIR/long_iso/passed/busco_passed.txt
  fi

else
  echo "ERROR: $OUT/final_annot.gff doesn't exist"
fi

}
# AGAT statistics
GFF_STATS(){
METHOD="braker3"

OUT=$ANNOT_DIR/$BASENAME_GENOME/${METHOD}_out

MODULES_LOAD AGAT
NAME="AGAT_STATS"

UNZ_GFFFASTA

if [ -f "$OUT/final_annot.gff" ] && [ -s "$OUT/final_annot.gff" ]
then

  if [ ! -f "$ANNOT_DIR/gff_stats/${BASENAME_GENOME}.txt" ] || [ ! -s "$ANNOT_DIR/gff_stats/${BASENAME_GENOME}.txt" ]
  then

    echo $'\n'"Extracting GFF stats at $(date)."

    agat_sp_statistics.pl \
    --gff $OUT/final_annot.gff \
    -o $ANNOT_DIR/gff_stats/${BASENAME_GENOME}.txt
    echo $'\n'"GFF stats obtained at $(date)."
  else
    echo "Output already exists."
  fi

else
  echo "ERROR: $OUT/final_annot.gff doesn't exist"
fi
}
#### Running commands
GENOME_ANNOT_BRAKER
FUNCT_ANNOT
GFF_STATS
BUSCO_P
TPM_CLEAN
