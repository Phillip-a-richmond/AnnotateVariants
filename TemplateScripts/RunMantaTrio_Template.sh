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
conda activate Manta


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
SIBLING_SAMPLEID=sibling_id_${GENOME}
PED=$FAMILY_ID.ped

PROBAND_BAM=${PROBAND_SAMPLEID}.dupremoved.sorted.bam
FATHER_BAM=${FATHER_SAMPLEID}.dupremoved.sorted.bam
MOTHER_BAM=${MOTHER_SAMPLEID}.dupremoved.sorted.bam

# MANTA stuff
ANALYSIS_DIR=$BAM_DIR/Manta

# ANNOTSV stuff
ANNOTSV=annotsv_dir
export ANNOTSV=annotsv_dir
GENELIST_BOOL=genelist_bool
GENELIST=genelist

## Step 1 - Make Config Script
configManta.py \
        --referenceFasta=$GENOME_FASTA \
        --runDir=$ANALYSIS_DIR \
        --bam $PROBAND_BAM \
        --bam $MOTHER_BAM \
        --bam $FATHER_BAM

## Step 2 - Execute config script
cd $ANALYSIS_DIR
./runWorkflow.py \
        -j $NSLOTS \
        -g 64

## Step 3 - Annotate
### With GeneList
if [ "$GENELIST_BOOL" = true ]; then
	$ANNOTSV/bin/AnnotSV -SVinputFile $ANALYSIS_DIR/results/variants/diploidSV.vcf.gz \
	        -genomeBuild $GENOME \
	        -overlap 80 \
	        -reciprocal yes \
	        -outputFile ${FAMILY_ID}-MANTA-annotsv-candidateGenes
fi
	
### Without GeneList
$ANNOTSV/bin/AnnotSV -SVinputFile $ANALYSIS_DIR/results/variants/diploidSV.vcf.gz \
        -genomeBuild $GENOME \
        -overlap 80 \
        -reciprocal yes \
        -outputFile ${FAMILY_ID}-MANTA-annotsv
