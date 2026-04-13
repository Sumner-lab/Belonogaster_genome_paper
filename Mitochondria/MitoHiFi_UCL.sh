#!/bin/bash
#SBATCH --job-name=MitoHiFi_UCL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
# "d-hh:mm:ss"
#SBATCH --mail-type=end
#SBATCH --mail-user=pk14885@bristol.ac.uk
#SBATCH --account=BISC020282

mkdir /user/work/pk14885/UCL_Genome_MitoHiFi
cd /user/work/pk14885/UCL_Genome_MitoHiFi

# ---------------------------------------------------------------------------------
# DOWNLOAD A REASONABLY CLOSE MITOCHONDRIAL GENOME TO ACT AS A REFERENCE:

# First, using MitoHiFi, we find a "closely related mitochondrial sequence" (https://github.com/marcelauliano/MitoHiFi).
## module load apptainer/1.3.1-ksax
## singularity run -B/user/work/pk14885/UCL_Genome_MitoHiFi/ /user/work/pk14885/mitohifi_master.sif findMitoReference.py --species "Polistes jokahamae" --outfolder /user/work/pk14885/UCL_Genome_MitoHiFi --min_length 14000
# This downloads the complete mitochondrial genome for a Polistes paper wasp (Polistes rothneyi), accession: https://www.ncbi.nlm.nih.gov/nuccore/NC_063693
# (Note that we use a Polistes wasp in the findMitoReference search above because Belonogaster leads to a failure to source a reasonably close genome.)
# Running that line downloads the P. rothneyi genome to the folder /user/work/pk14885/UCL_Genome_MitoHiFi/ (i.e., you should see the files LC519884.1.fasta  LC519884.1.gb)

# ---------------------------------------------------------------------------------
# RUN MitoHiFi ON BLUEPEBBLE:

module load apptainer/1.3.1-ksax

singularity run -B/user/work/pk14885/UCL_Genome_MitoHiFi/ /user/work/pk14885/mitohifi_master.sif mitohifi.py \
        -r /user/work/pk14885/UCL_Genome_MitoHiFi/DTG-DNA-1074.r64296e222347B01.subreads_ccs.fastq.mini.gz \
        -f LC519884.1.fasta \
        -g LC519884.1.gb -t 8 -o 5 

 # ---------------------------------------------------------------------------------

# The result of this script is the annotated mitochondrial genome for the focal wasp.

# ---------------------------------------------------------------------------------
# CALCULATING MEAN DEPTH:

# To calculate mean depth, we use minimap2:
## minimap2 -ax map-hifi final_mitogenome.fasta DTG-DNA-1074.r64296e222347B01.subreads_ccs.fastq.gz | samtools sort -o mito.coverage.bam
## samtools index mito.coverage.bam
## samtools flagstat mito.coverage.bam
## samtools depth mito.coverage.bam | awk '{sum+=$3} END {print sum/NR}'
















