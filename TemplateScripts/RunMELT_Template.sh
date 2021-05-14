#!/bin/bash

#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########

# Load software environment for MELT (bowtie2/samtools)
MINICONDA3_DIR=miniconda3_dir
source $MINICONDA3_DIR/opt/miniconda3/etc/profile.d/conda.sh
conda activate MELT

# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR
WORKING_DIR=working_dir

## Set working space
mkdir -p $WORKING_DIR
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=genome_build
FASTA_DIR=fasta_dir
FASTA_FILE=fasta_file
GENOME_FASTA=$FASTA_DIR/$FASTA_FILE

BAM_DIR=$WORKING_DIR
FAMILY_ID=family_id
PROBAND_SAMPLEID=proband_id_${GENOME}
MOTHER_SAMPLEID=mother_id_${GENOME}
FATHER_SAMPLEID=father_id_${GENOME}
PED=$FAMILY_ID.ped

PROBAND_BAM=${PROBAND_SAMPLEID}.dupremoved.sorted.bam
FATHER_BAM=${FATHER_SAMPLEID}.dupremoved.sorted.bam
MOTHER_BAM=${MOTHER_SAMPLEID}.dupremoved.sorted.bam

ANALYSIS_DIR=${WORKING_DIR}MELT/
MELT_DIR=melt_dir
MEI_LIST=${MELT_DIR}/me_refs/GRCh38/mei_list.txt
GENE_ANNO=$MELT_DIR/add_bed_files/GRCh38/GRCh38.genes.bed
BAM1=$PROBAND_BAM
BAM2=$MOTHER_BAM
BAM3=$FATHER_BAM

## Step 1 - Preprocess
java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
        -bamfile $BAM1 \
	-h $GENOME_FASTA

java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
        -bamfile $BAM2 \
	-h $GENOME_FASTA

java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
        -bamfile $BAM3 \
	-h $GENOME_FASTA

# Step 2 - Individual Analysis
java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
	-w $ANALYSIS_DIR \
	-t $MEI_LIST \
	-c 30 \
	-bamfile $BAM1 \
	-h $GENOME_FASTA

java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
        -w $ANALYSIS_DIR \
        -t $MEI_LIST \
        -c 30 \
        -bamfile $BAM2 \
        -h $GENOME_FASTA

java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
        -w $ANALYSIS_DIR \
        -t $MEI_LIST \
        -c 30 \
        -bamfile $BAM3 \
        -h $GENOME_FASTA

# Step 3 Group Discovery
java -Xmx4G -jar ${MELT_DIR}MELT.jar GroupAnalysis \
	-discoverydir $ANALYSIS_DIR \
	-w $ANALYSIS_DIR \
	-t $MEI_LIST \
	-h $GENOME_FASTA \
	-n $GENE_ANNO 
	

# Step 4 Genotype
java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
	-bamfile $BAM1 \
	-t $MEI_LIST \
	-h $GENOME_FASTA \
	-w $ANALYSIS_DIR \
	-p $ANALYSIS_DIR 

java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
        -bamfile $BAM2 \
        -t $MEI_LIST \
        -h $GENOME_FASTA \
        -w $ANALYSIS_DIR \
        -p $ANALYSIS_DIR 

java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
        -bamfile $BAM3 \
        -t $MEI_LIST \
        -h $GENOME_FASTA \
        -w $ANALYSIS_DIR \
        -p $ANALYSIS_DIR 


# Make VCF
java -Xmx2G -jar ${MELT_DIR}MELT.jar MakeVCF \
	-genotypingdir $ANALYSIS_DIR \
	-h $GENOME_FASTA \
	-t $MEI_LIST \
	-w $ANALYSIS_DIR \
	-p $ANALYSIS_DIR \
	-o $ANALYSIS_DIR


