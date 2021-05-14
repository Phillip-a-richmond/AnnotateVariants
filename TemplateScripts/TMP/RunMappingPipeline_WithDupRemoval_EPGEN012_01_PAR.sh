#!/bin/bash

#SBATCH --partition=defq

#SBATCH --mail-user=malavika.hebbar@cw.bc.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=160G
# Try to keep RAM=cores*4, since there is roughly 4G RAM/core

#SBATCH --cpus-per-task=40
#SBATCH --time=24:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error



##########
# Set up #
##########

## Get the tools we need, from a conda environment within WASSERMAN_SOFTWARE
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR

# Set working space
WORKING_DIR=/mnt/scratch/Precision/EPGEN/PROCESS/EPGEN012_PAR/
mkdir -p $WORKING_DIR
RAW_DIR=/mnt/scratch/Precision/EPGEN/RAW/EPGEN012/
NSLOTS=$SLURM_CPUS_PER_TASK

SAMPLE=EPGEN012-01

echo $SAMPLE
echo "${SAMPLE}_1.fastq.gz"
echo "${SAMPLE}_2.fastq.gz"

FASTQR1=$RAW_DIR${SAMPLE}_1.fastq.gz
FASTQR2=$RAW_DIR${SAMPLE}_2.fastq.gz

ls $FASTQR1
ls $FASTQR2


# STEP 1: MAP AGAINST GRCH38 #

# Genome Information
#### GRCh38 ####
BWA_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa
# Set SAMPLE_ID to have GRCh38
SAMPLE_ID=${SAMPLE}_GRCh38

##################
# Initialize Log #
##################
# Timing from here to end of script
FullRunStart=`date +%s`

LOGFILE=${SAMPLE_ID}_logfile.csv
rm ${SAMPLE_ID}_logfile.csv
touch ${SAMPLE_ID}_logfile.csv
# This logfile will be written to in csv format.
# Columns for log file:
echo "SampleID,Operation,Runtime">> $LOGFILE

#############
# Map reads #
#############
cd $WORKING_DIR

Start=`date +%s`

bwa mem $BWA_INDEX \
        -t $NSLOTS \
        -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" \
	-Y -K 100000000 \
        $FASTQR1 \
        $FASTQR2 \
        > $WORKING_DIR$SAMPLE_ID.sam

End=`date +%s`
runtime=$((End-Start))
echo "$SAMPLE_ID,BWAmem,$runtime" >> $LOGFILE
echo "BWA mem ran in $runtime"


############
# Samtools #
############
#convert to binary and index
Start=`date +%s`

samtools view -@ $NSLOTS -ubS $WORKING_DIR${SAMPLE_ID}.sam \
        | samtools sort - -@ $NSLOTS  -T $WORKING_DIR${SAMPLE_ID}.sorted -O BAM -o $WORKING_DIR${SAMPLE_ID}.sorted.bam
samtools index $WORKING_DIR${SAMPLE_ID}.sorted.bam

End=`date +%s`
runtime=$((End-Start))
echo "$SAMPLE_ID,Samtools,$runtime" >> $LOGFILE
echo "Samtools ran in $runtime"

rm $WORKING_DIR${SAMPLE_ID}.sam

######################
# Picard Dup Removal #
######################

# Unload bcchr, and load cvmfs
# unload_bcchr
source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
# load cvmfs
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load picard

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
	R=$BWA_INDEX \
	I=$WORKING_DIR${SAMPLE_ID}.sorted.bam \
	O=$WORKING_DIR${SAMPLE_ID}.dupremoved.sorted.bam \
	REMOVE_DUPLICATES=false \
	M=$WORKING_DIR${SAMPLE_ID}.duplicateMetrics.txt

## Get the tools we need, from a conda environment within WASSERMAN_SOFTWARE
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

samtools index $WORKING_DIR${SAMPLE_ID}.dupremoved.sorted.bam

exit
