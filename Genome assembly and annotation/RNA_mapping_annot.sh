#!/bin/bash

########### Do not modify the following
#-------------------------------------------------------------------------------------------------------------

####################################################
# Error handler
####################################################
EH() {
# exit when any command fails
#set -e
# keep track of the last executed command
#trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
#trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT
echo ""
}
####################################################
# HELP
####################################################
HELP() {
  echo ""
  echo "Hello, this program will execute 2 pipelines for RNA-seq data cleaning and mapping for genome/transcriptome annotation."
  echo ""
  echo "First, you have to choose between pipeline 1 (trimming of the raw data) and pipeline 2 (RNAseq mapping on ref genome)."
  echo "Each pipeline is divided into steps:"
  echo "Pipeline 1: step 1 (trimmomatic), step 2 (quality control check) and all (steps1 and 2)."
  echo "Pipeline 2: step 1 (mapping), step 2 (bam sorting) and all (steps1 and 2)."
  echo ""
  echo "syntax = 01_RNA_mappinf_annot.sh -c config/file/path -p 1 -s all"
  echo "    -c      Path to your config file."
  echo "    -p      Pipeline to execute (1 or 2)."
  echo "    -s      Step in the pipeline to execute (values will be 1,2 or all for -p 1 ; 1,2 or all for -p 2)."
  echo "    -h      Will display this help."
  echo "The syntax above will perform pipeline 1 (trimming and filtering) with all its steps."
  echo ""
  echo "If you have any problems, you can contact me at this address alexandra.cdea@gmail.com"
  echo ""
}

####################################################
# Option handler
####################################################
while getopts ":hc:p:s:" option; do
   case $option in
      h) # display Help
         HELP
         exit;;
      c) # path
         CONFIG=${OPTARG};;
      p) # pipeline 
         PIPE=${OPTARG};;
      s) # step
         STEP=${OPTARG};;
     \?) # Invalid option
         echo "Error: Invalid option"
         echo "Use the syntax: 01_RNA_mapping_annot.sh -c config/file/path -p 1 -s all"
         echo "use -h option to print the help page"
         exit 2;;
   esac
done
####################################################

if [ ! -f "$CONFIG" ]
then
echo "Can't access $CONFIG config file! Please, verify your path."
HELP
echo "Now exiting..."
exit 1
fi

if [ -z "$PIPE" ] || [ -z "$STEP" ]
then
HELP
exit 2
fi

## working directory

WORK_DIR=$(grep "WORK_DIR:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")  # working dir path

## Defaults

TEMP_DIR=$WORK_DIR/temp
RAW_DATA_extension=.fastq.gz
RAW_DATA_forward_extension=_1.fastq.gz
RAW_DATA_reverse_extension=_2.fastq.gz
LOG_FILE=$WORK_DIR/log
GENOME_EXTENSION=.fna
N_THREADS=4
TRIM_ADAPTORS_PE=TruSeq3-PE.fa
TRIM_ADAPTORS_SE=TruSeq3-SE.fa

##
TRIM_ADAPTORS_PATH=$(grep "TRIM_ADAPTORS_PATH:" $CONFIG | grep -oP "(?<=: ).*(?=#)" | sed "s/^ *//;s/ *$//") # path to the trimmomatic adaptor directory
##

## Defaults
TRIMMOMATIC_OPTIONS=$(echo "-threads $N_THREADS -phred33 ILLUMINACLIP:$TRIM_ADAPTORS_PATH/${TRIM_ADAPTORS_PE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
HISAT_OPTIONS=$(echo "--phred33 --rna-strandness RF --dta -t")

## Directories

TEMP_DIR=$(grep "TEMP_DIR:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # path to temporary directory
RAW_DATA_DIR=$(grep "RAW_DATA_DIR:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # raw reads path
RAW_DATA_extension=$(grep "RAW_DATA_extension:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # extension / file type of your raw data
RAW_DATA_forward_extension=$(grep "RAW_DATA_forward_extension:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # forward file extension of your RNA-seq files
RAW_DATA_reverse_extension=$(grep "RAW_DATA_reverse_extension:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # reverse file extension of your RNA-seq files
LOG_FILE=$(grep "LOG_FILE:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # path to your log directory

GENOME_FILE=$(grep "GENOME_FILE:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # Path to your genome file (used in pipeline2)
GFF_FILE=$(grep "GFF_FILE:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")
BASENAME=$(grep "BASENAME:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")
R_SCRIPT=$(grep "SCRIPT_PATH:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")
GENOME_EXTENSION=$(grep "GENOME_EXTENSION:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # Extension of you genome file
TO_MERGE=( $(grep "TO_MERGE:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") ) # If you want to merge some files (put the commun extension of files to merge, ex: P1_S1_fwd.fq P1_S2_fwd.fq, P1 would be precised if you want to merge by P1 samples)

TRIMMOMATIC_PATH=$(grep "TRIMMOMATIC_PATH:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")

## options

N_THREADS=$(grep "N_THREADS:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # number of threads you want to use
PAIRED_END=$(grep "PAIRED_END:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # Are the raw reads paired-end? yes/no
TRIM_ADAPTORS_PE=$(grep "TRIM_ADAPTORS_PE:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") # Adaptor file you want to use. PE for paired-end, SE for single-end
TRIM_ADAPTORS_SE=$(grep "TRIM_ADAPTORS_SE:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")
TRIMMOMATIC_OPTIONS=$(echo $(grep "TRIMMOMATIC_OPTIONS:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") | sed "s+\$N_THREADS+$N_THREADS+1;s+\$TRIM_ADAPTORS_PATH+$TRIM_ADAPTORS_PATH+;s/${TRIM_ADAPTORS_SE}/${TRIM_ADAPTORS_SE}/;s+\${TRIM_ADAPTORS_PE}+${TRIM_ADAPTORS_PE}+" ) # trimmomatic options, do not precise PE/SE, -trimolog and -summary options

HISAT_OPTIONS=$(echo $(grep "HISAT_OPTIONS:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")) # hisat options, do not precise -1, -2, -S,and --novel-splicesite-outfile options
#JOB_OPTIONS=$(echo $(grep "JOB_OPTIONS:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//") | sed "s+\$N_THREADS+$N_THREADS+1" )

MODULES=($(grep "MODULES:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//"))
MODULES_STEP4=($(grep "MODULES_STEP4:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//"))
R_ENV=$(grep "R_ENV:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")
ENV_ACTIVATION=$(grep "ENV_ACTIVATION:" $CONFIG | grep -oP "(?<=:).*(?=#)" | sed "s/^ *//;s/ *$//")

echo "Starting transcriptome alignment pipeline at $(date)"
echo

if [ -z "$WORK_DIR" ] || [ -z "$RAW_DATA_DIR" ] || [ -z "$TRIM_ADAPTORS_PATH" ] || [ -z "$TRIMMOMATIC_PATH" ]
then 
echo "Error: something is wrong with your config file! Please, use a config file with the right format."
echo "You must precise at least the working directory path, the raw data directory path, the genome path and the path to the adapter sequences."
echo
echo "Mandatory fields:"
echo "WORK_DIR"
echo "RAW_DATA_DIR"
echo "TRIM_ADAPTORS_PATH"
echo "TRIMMOMATIC_PATH"
echo
echo "Exiting..."
exit 1
fi

echo "Options: -c $CONFIG -p $PIPE -s $STEP"
echo "Working directory: $WORK_DIR"
echo "Raw data directory: $RAW_DATA_DIR"
echo "Genome path: $GENOME_FILE"
echo "Adaptor path: $TRIM_ADAPTORS_PATH"
echo "Temprary directory: $TEMP_DIR"
echo "Raw data exetnsion: \"$RAW_DATA_extension\" | forward: \"$RAW_DATA_forward_extension\" | reverse: \"$RAW_DATA_reverse_extension\""
echo "Log file pat: $LOG_FILE"
echo "Genome extension: \"$GENOME_EXTENSION\""
echo "Number of threads: $N_THREADS"
echo "Adaptors to concider: $TRIM_ADAPTORS_PE $TRIM_ADAPTORS_SE"
echo "Path to trimmomatic: $TRIMMOMATIC_PATH"
echo "Trimmomatic options: $TRIMMOMATIC_OPTIONS"
echo "Merging: $TO_MERGE"
echo "job submission management: $WORK_MANAGEMENT_NAME"
echo "job options: $JOB_OPTIONS"
echo

if [ ! -f "$TRIMMOMATIC_PATH" ]
then
  echo "Warning: $TRIMMOMATIC_PATH is not a file or does not exist."
fi

#-------------------------------------------------------------------------------------------------------------
####################################################
## Pipeline 1
####################################################

START1() {
if [ ! -d "$WORK_DIR/01_trimmed_filtered_data" ]; then
  echo "Creating output direrctory 01_trimmed_filtered_data ."
  mkdir -p $WORK_DIR/01_trimmed_filtered_data
fi

OUT_DIR=$WORK_DIR/01_trimmed_filtered_data

}

## trimmomatic on raw reads
STEP1() {

echo "Starting Trimmomatic."

START1

if [ ! -d "$OUT_DIR/01_trimmomatic" ]; then
  echo "Creating 01_trimmed_data directory in $OUT_DIR"
  mkdir -p $OUT_DIR/01_trimmomatic
fi

TRIM_OUT_DIR=$OUT_DIR/01_trimmomatic

if [ ! -f "$TRIMMOMATIC_PATH" ]
then
echo "Error: $TRIMMOMATIC_PATH is not a file or does not exist."
echo "Exiting..."
exit 1
fi

for fastq_files in $(ls $RAW_DATA_DIR/*${BASENAME}*${RAW_DATA_extension} | xargs -n1 basename)
do
  if [ ! -z "$(echo $fastq_files | grep "${RAW_DATA_forward_extension}$")" ]
  then
    echo "Filtering on ${fastq_files} and ${fastq_files/$RAW_DATA_forward_extension/$RAW_DATA_reverse_extension} at $(date)."
    java -jar $TRIMMOMATIC_PATH PE -trimlog $TRIM_OUT_DIR/${fastq_files/$RAW_DATA_forward_extension/}.log \
    -summary $TRIM_OUT_DIR/${fastq_files/$RAW_DATA_forward_extension/}.stats.txt \
    $RAW_DATA_DIR/${fastq_files} \
    $RAW_DATA_DIR/${fastq_files/$RAW_DATA_forward_extension/$RAW_DATA_reverse_extension} \
    $TRIM_OUT_DIR/paired.${fastq_files} \
    $TRIM_OUT_DIR/unpaired.${fastq_files} \
    $TRIM_OUT_DIR/paired.${fastq_files/$RAW_DATA_forward_extension/$RAW_DATA_reverse_extension} \
    $TRIM_OUT_DIR/unpaired.${fastq_files/$RAW_DATA_forward_extension/$RAW_DATA_reverse_extension} \
    $TRIMMOMATIC_OPTIONS
    echo "Filtering finished for ${fastq_files} and ${fastq_files/$RAW_DATA_forward_extension/$RAW_DATA_reverse_extension}."
  elif [ ! -z "$(echo $fastq_files | grep "${RAW_DATA_reverse_extension}$")" ]
  then
    echo "reverse read file"
  else
    echo "Filtering on ${fastq_files} at $(date)."
    java -jar $TRIMMOMATIC_PATH SE -trimlog $TRIM_OUT_DIR/${fastq_files/$RAW_DATA_extension/}.log \
    -summary $TRIM_OUT_DIR/${fastq_files/$RAW_DATA_extension/}.stats.txt \
    $RAW_DATA_DIR/${fastq_files} \
    $TRIM_OUT_DIR/trimmed.${fastq_files} \
    ${TRIMMOMATIC_OPTIONS/$TRIM_ADAPTORS_PE/$TRIM_ADAPTORS_SE}
    echo "Filtering finished for ${fastq_files}."
  fi
done

echo "Trimmomatic finished at $(date)."
}

## Quality control
STEP2() {
echo "Starting MultiQC."

START1

if [ ! -d "$OUT_DIR/02_quality_control" ] ; then
  echo "Creating 02_quality_control directory."
  mkdir -p $OUT_DIR/02_quality_control
fi

multiQC_OUT_DIR=$OUT_DIR/02_quality_control

IN_DIR=$OUT_DIR/01_trimmomatic

if [ ! -x "$(command -v fastqc)" ] 
then
echo "Error: fastqc is not in your PATH."
echo "Exiting..."
exit 126
fi

if [ ! -x "$(command -v multiqc)" ]
then
echo "Error: multiqc is not in your PATH."
echo "Exiting..."
exit 126
fi

for trimmomatic_seq in $(ls $IN_DIR/paired.*${BASENAME}*${RAW_DATA_extension} | grep -v "unpaired" | xargs -n1 basename)
do
  dir=${trimmomatic_seq/${RAW_DATA_extension}/}
  mkdir -p $multiQC_OUT_DIR/$dir/fastQC_report/
  fastqc -t $N_THREADS -o $multiQC_OUT_DIR/$dir/fastQC_report $IN_DIR/$trimmomatic_seq
done
multiqc $multiQC_OUT_DIR/*/fastQC_report/ --force --outdir $multiQC_OUT_DIR

echo "MultiQC finished at $(date)."
}

#-------------------------------------------------------------------------------------------------------------
####################################################
## Pipeline 2
####################################################

START2() {
if [ ! -d "$WORK_DIR/02_mapping_assembly" ]; then
  echo "Creating output direrctory 02_mapping_assembly ."
  mkdir -p $WORK_DIR/02_mapping_assembly
fi

OUT_DIR=$WORK_DIR/02_mapping_assembly
}

## Hisat : mapping RNA reads on genome
STEP21() {
echo "Starting Hisat2."

START2

IN_DIR=$WORK_DIR/01_trimmed_filtered_data

if [ ! -d "$OUT_DIR/01_hisat_mapping" ] ; then
  echo "Creating 01_hisat_mapping."
  mkdir -p $OUT_DIR/01_hisat_mapping
fi

HISAT_OUT_DIR=$OUT_DIR/01_hisat_mapping

# sam outdir
if [ ! -d "$OUT_DIR/02_bam_files" ] ; then
  echo "Creating 02_bam_files."
  mkdir -p $OUT_DIR/02_bam_files
fi

SAMTOOLS_OUT_DIR=$OUT_DIR/02_bam_files

if [ ! -x "$(command -v hisat2)" ]
then
echo "Error: hisat2 is not in your PATH."
echo "Exiting..."
exit 126
fi

if [ ! -x "$(command -v samtools)" ]
then
echo "Error: samtools is not in your PATH."
echo "Exiting..."
exit 126
fi

echo "Building index for Hisat."

if [ ! -d "$HISAT_OUT_DIR/index" ] ; then
  echo "Creating index directory."
  mkdir -p $HISAT_OUT_DIR/index
fi

if [ -z "$(ls $HISAT_OUT_DIR/index/ | grep $(basename ${GENOME_FILE/$GENOME_EXTENSION/}))" ]
then
  echo "Building index..."
  hisat2-build -p 7 $GENOME_FILE $HISAT_OUT_DIR/index/$(basename ${GENOME_FILE/$GENOME_EXTENSION/})
fi
INDEX_PATH=$HISAT_OUT_DIR/index/$(basename ${GENOME_FILE/$GENOME_EXTENSION/})

echo "Index built at $(date)."

for trimmed_seq in $(ls $IN_DIR/01_trimmomatic/*${BASENAME}*${RAW_DATA_extension} | grep -v "unpaired" | xargs -n1 basename); do
  if [ ! -z "$(echo $trimmed_seq | grep "paired" | grep "${RAW_DATA_forward_extension}$")" ]
  then
    echo "Starting Hisat mapping of ${trimmed_seq} on $(basename ${GENOME_FILE/$GENOME_EXTENSION/}) genome."
    hisat2 -x $INDEX_PATH \
    -1 $IN_DIR/01_trimmomatic/${trimmed_seq} \
    -2 $IN_DIR/01_trimmomatic/${trimmed_seq/${RAW_DATA_forward_extension}/${RAW_DATA_reverse_extension}} \
    -S $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_forward_extension}/}-Hisat.sam \
    --novel-splicesite-outfile $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_forward_extension}/}-Hisat.junctions \
    $HISAT_OPTIONS
## SAM to BAM
    echo "Converting $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_forward_extension}/}-Hisat.sam into bam at $(date)."
    samtools view -S -bo $SAMTOOLS_OUT_DIR/${trimmed_seq/${RAW_DATA_forward_extension}/}-Hisat.bam \
    $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_forward_extension}/}-Hisat.sam
    if [ -f "$SAMTOOLS_OUT_DIR/${trimmed_seq/${RAW_DATA_forward_extension}/}-Hisat.bam" ] && [ $(stat -c %s "$SAMTOOLS_OUT_DIR/${trimmed_seq/${RAW_DATA_forward_extension}/}-Hisat.bam") -gt 100000 ]
    then
      rm $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_forward_extension}/}-Hisat.sam
    fi
  elif [ ! -z "$(echo $trimmed_seq | grep "paired" | grep "${RAW_DATA_reverse_extension}$")" ]
  then
    echo "reverse read file"
  else
    echo "Starting Hisat mapping of $(trimmed_seq) on $(basename ${GENOME_FILE/$GENOME_EXTENSION/}) genome."
    hisat2 -x $INDEX_PATH \
    -U $IN_DIR/01_trimmomatic/${trimmed_seq} \
    -S $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_extension}/}-Hisat.sam \
    --novel-splicesite-outfile $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_extension}/}-Hisat.junctions \
    $HISAT_OPTIONS
## SAM to BAM
    echo "Converting $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_extension}/}-Hisat.sam into bam at $(date)."
    samtools view -S -bo $SAMTOOLS_OUT_DIR/${trimmed_seq/${RAW_DATA_extension}/}-Hisat.bam \
    $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_extension}/}-Hisat.sam
    if [ -f "$SAMTOOLS_OUT_DIR/${trimmed_seq/${RAW_DATA_extension}/}-Hisat.bam" ] && [ $(stat -c %s "$SAMTOOLS_OUT_DIR/${trimmed_seq/${RAW_DATA_extension}/}-Hisat.bam") -gt 100000 ]
    then
      rm $HISAT_OUT_DIR/${trimmed_seq/${RAW_DATA_extension}/}-Hisat.sam
    fi
  fi
done

echo "Hisat finished at $(date)."
}

## BAM sorting

STEP22() {
echo "Starting bam sorting."

START2

if [ ! -d "$OUT_DIR/03_sorted_bam_files" ] ; then
  echo "Creating 03_sorted_bam_files."
  mkdir -p $OUT_DIR/03_sorted_bam_files
fi

SAMTOOLS_IN_DIR=$OUT_DIR/02_bam_files
SAMTOOLS_OUT_DIR=$OUT_DIR/03_sorted_bam_files

if [ ! -x "$(command -v samtools)" ]
then
echo "Error: samtools is not in your PATH."
echo "Exiting..."
exit 126
fi


for bam_files in $(ls $SAMTOOLS_IN_DIR/*${BASENAME}*.bam | xargs -n1 basename); do
  echo "Sorting $bam_files at $(date)."
  samtools sort -T $TEMP_DIR \
  -o $SAMTOOLS_OUT_DIR/${bam_files/.bam/_sorted.bam} \
  $SAMTOOLS_IN_DIR/$bam_files
  samtools index $SAMTOOLS_OUT_DIR/${bam_files/.bam/_sorted.bam}
done
echo "Bam sorting finished at $(date)."
}

#-------------------------------------------------------------------------------------------------------------

if [ "$PIPE" == 1 ] && [ "$STEP" == 1 ]; then
  {
  EH
  STEP1
  }
# >& $LOG_FILE/01_trimmomatic.log
elif [ "$PIPE" == 1 ] && [ "$STEP" == 2 ]; then
  {
  EH
  STEP2 
  }
# >& $LOG_FILE/02_MultiQC.log
elif [ "$PIPE" == 1 ] && [ "$STEP" == "all" ]; then
  {
  EH
  STEP1
  STEP2
  }
# >& $LOG_FILE/all_pipeline-1.log
elif [ "$PIPE" == 2 ] && [ "$STEP" == 1 ]; then
  {
  EH
  STEP21
  }
# >& $LOG_FILE/02-1_hisat.log
elif [ "$PIPE" == 2 ] && [ "$STEP" == 2 ]; then
  {
  EH
  STEP22 
  }
# >& $LOG_FILE/02-2_bam_sorting.log
elif [ "$PIPE" == 2 ] && [ "$STEP" == "all" ]; then
  {
  EH
  STEP21
  STEP22
  } 
# >& $LOG_FILE/all_pipeline-2.log
else
  HELP
fi
