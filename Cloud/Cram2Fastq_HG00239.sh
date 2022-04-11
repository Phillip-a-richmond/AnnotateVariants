#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=160G
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########

# Load environment (need this for outlier with ehdn)
source /mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/opt/miniconda3/etc/profile.d/conda.sh 
conda activate /mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/opt/AnnotateVariantsEnvironment/


# Number of threads
NSLOTS=40

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR
WORKING_DIR=/mnt/common/OPEN_DATA/ThousandGenomes/

## Set working space
mkdir -p $WORKING_DIR
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=GRCh38
FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/
FASTA_FILE=GRCh38_full_analysis_set_plus_decoy_hla.fa

INCRAM=HG00239.final.cram
INBAM=HG00239.final.bam
INBAM_SORTED=HG00239.namesorted.bam
OUTFQ1=HG00239.R1.fastq
OUTFQ2=HG00239.R2.fastq


# Samtools view
samtools view -hbu -@ $NSLOTS -T ${FASTA_DIR}/$FASTA_FILE $INCRAM | \
	samtools sort -n -@ $NSLOTS -o $INBAM_SORTED -

# Convert to fastq
samtools fastq \
	-@ 78 \
	-1 $OUTFQ1 -2 $OUTFQ2 \
	-0 /dev/null \
	-s /dev/null \
	-n $INBAM_SORTED

gzip $OUTFQ1
gzip $OUTFQ2

