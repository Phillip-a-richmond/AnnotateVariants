#!/bin/bash

#SBATCH --partition=defq

#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=80G
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

## Job Array stuff
#  #SBATCH --array=0-30%3

##########
# Set up #
##########

## Get the tools we need, from a conda environment within WASSERMAN_SOFTWARE
source /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda activate STAR

# Threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Set working space
WORKING_DIR=/mnt/common/OPEN_DATA/BuildControlCohort/RNA/
mkdir -p $WORKING_DIR
RAW_DIR=/mnt/common/OPEN_DATA/BuildControlCohort/RNA/

# Define sample id from set of fastq files, based on the job array index
SLURM_ARRAY_TASK_ID=0
Files=(${RAW_DIR}*_1.fastq.gz)
IFS='/' read -a array <<< ${Files[$SLURM_ARRAY_TASK_ID]}
SampleR1Fastq=${array[-1]}
IFS='_' read -a array2 <<< "${SampleR1Fastq}"
SAMPLE=${array2[0]}_mRNA

echo $SAMPLE
echo "${SAMPLE}_1.fastq.gz"
echo "${SAMPLE}_2.fastq.gz"

FASTQR1=$RAW_DIR${SAMPLE}_1.fastq.gz
FASTQR2=$RAW_DIR${SAMPLE}_2.fastq.gz

ls $FASTQR1
ls $FASTQR2
exit
# STEP 1: MAP AGAINST GRCH38 #

# Genome Information
#### GRCh38 ####
STAR_INDEX=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GENCODE/STAR/GRCh38_GENCODEv38_TranscriptAware
STAR_SPLICESITES=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GENCODE/STAR/gencode.v38.annotation.gtf.hisat2.splicesites.tsv

# Set SAMPLE_ID to have GRCh38, STAR
SAMPLE_ID=${SAMPLE}_STAR_GRCh38


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
hisat2 \
	-p $NSLOTS \
	--known-splicesite-infile $STAR_SPLICESITES \
	-x $STAR_INDEX \
	-1 $FASTQR1 \
	-2 $FASTQR2 \
	-S $SAMPLE_ID.sam


End=`date +%s`
runtime=$((End-Start))
echo "$SAMPLE_ID,STAR,$runtime" >> $LOGFILE
echo "STAR ran in $runtime"


############
# Samtools #
############
#convert to binary and index
Start=`date +%s`

samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
        | samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

End=`date +%s`
runtime=$((End-Start))
echo "$SAMPLE_ID,Samtools,$runtime" >> $LOGFILE
echo "Samtools ran in $runtime"

rm $WORKING_DIR$SAMPLE_ID'.sam'


# 

exit



