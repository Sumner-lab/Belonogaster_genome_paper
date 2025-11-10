#!/bin/bash
#SBATCH --job-name=MitoZ_UCL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=0-06:00:00
# "d-hh:mm:ss"
#SBATCH --mail-type=end
#SBATCH --mail-user=pk14885@bristol.ac.uk
#SBATCH --account=BISC020282

cd /user/work/pk14885/UCL_Genome_MitoZ

# ---------------------------------------------------------------------------------
# RUN MITOZ ON BLUEPEBBLE (OPTION 1):

# Running it through Singularity on BluePebble (Bristol's HPC):

module load apptainer/1.3.1-ksax

singularity run -B/user/work/pk14885/UCL_Genome_MitoZ/ /user/work/pk14885/MitoZ_v3.6.sif mitoz all \
           --fq1 /user/work/pk14885/UCL_Genome_MitoZ/DTG-DNA-1074.r64296e222347B01.subreads_ccs.fastq.mini.gz \
           --outprefix UCL_Genome_ \
           --genetic_code 5 \
           --clade Arthropoda \
           --fastq_read_length 150 \
           --data_size_for_mt_assembly 0.5 \
           --assembler megahit \
           --kmers_megahit 21 29 39 59 79 99 \
           --thread_number 1 \
           --memory 30 \
           --requiring_taxa Arthropoda

 # ---------------------------------------------------------------------------------

# The result of this script is the annotated mitochondrial genome for the focal wasp.

# ---------------------------------------------------------------------------------
