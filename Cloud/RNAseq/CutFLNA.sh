#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=80G
#SBATCH --cpus-per-task=10
#SBATCH --time=86:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

## Job Array stuff
#SBATCH --array=1-28%6

##########
# Set up #
##########

# Load singularity
module load singularity
BIN_VERSION="1.1.0"

# Load env for bcftools
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment


# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Used for debuggin
#SLURM_ARRAY_TASK_ID=0

# Set working space
WORKING_DIR=/mnt/common/OPEN_DATA/BuildControlCohort/WGS/
mkdir -p $WORKING_DIR
BAM_DIR=/mnt/common/OPEN_DATA/BuildControlCohort/WGS/

# Define sample id from set of fastq files, based on the job array index
Files=(${BAM_DIR}*cram)
IFS='/' read -a array <<< ${Files[$SLURM_ARRAY_TASK_ID]}
SampleBam=${array[-1]}
IFS='.' read -a array2 <<< "${SampleBam}"
SAMPLE=${array2[0]}

echo $SAMPLE
echo "${SAMPLE}.final.cram"

BAM=$BAM_DIR/$SAMPLE.final.cram

ls $BAM

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=GRCh38
FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/
FASTA_FILE=GRCh38_full_analysis_set_plus_decoy_hla.fa 

SEQ_TYPE=WGS
BAM_DIR=$WORKING_DIR
TMP_DIR=$WORKING_DIR/tmp
mkdir -p $TMP_DIR

# FLNA bed file
TARGET=/mnt/common/OPEN_DATA/BuildControlCohort/FLNA.bed
# Index CRAM
samtools view \
	--write-index \
	$BAM \
	-@ $NSLOTS \
	-T $FASTA_DIR/$FASTA_FILE \
	-C \
	-L $TARGET \
	-o ${SAMPLE}_FLNA-target-region.cram

