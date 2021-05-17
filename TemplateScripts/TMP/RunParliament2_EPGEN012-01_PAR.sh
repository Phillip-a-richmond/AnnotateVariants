#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=160G
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########

# Load singularity
module load singularity

# Load env for bcftools
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Pull latest version, if you already have it, this will be skipped
#singularity pull docker://dnanexus/parliament2:v0.1.11-0-gb492db6d
# This SIF is here:
PARLIAMENT2_SIF=/mnt/common/Precision/Parliament2/parliament2_v0.1.11-0-gb492db6d.sif

# Arrange your working space
SAMPLE=EPGEN012-02
INBAM=${SAMPLE}_GRCh38.dupremoved.sorted.bam
WORKING_DIR=/mnt/scratch/Precision/EPGEN/PROCESS/EPGEN012_PAR/
INPUT_DIR=${WORKING_DIR}/Parliament2-${SAMPLE}/Input
OUTPUT_DIR=${WORKING_DIR}/Parliament2-${SAMPLE}/Output
TMP_DIR=$INPUT_DIR/tmp

mkdir -p ${WORKING_DIR}/Parliament2-${SAMPLE}/
mkdir -p $OUTPUT_DIR
mkdir -p $INPUT_DIR
mkdir -p $TMP_DIR

FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/
FASTA_FILE=GRCh38_full_analysis_set_plus_decoy_hla.fa

# Copy files to your input dir
# Haven't figured out a better way for this yet with these damn containers
if [ ! -f $INPUT_DIR/$FASTA_FILE ]; then
	cp $FASTA_DIR/$FASTA_FILE $INPUT_DIR/$FASTA_FILE
	cp $FASTA_DIR/$FASTA_FILE.fai $INPUT_DIR/$FASTA_FILE.fai
fi

if [ ! -f $INPUT_DIR/$INBAM ]; then
	cp $WORKING_DIR/$INBAM $INPUT_DIR
	cp $WORKING_DIR/$INBAM.bai $INPUT_DIR
fi

cd $INPUT_DIR
LANG= singularity run --workdir $TMP_DIR \
        -B ${INPUT_DIR}:/home/dnanexus/in:rw \
        -B ${OUTPUT_DIR}:/home/dnanexus/out:rw \
	$PARLIAMENT2_SIF \
	--bam $INBAM \
	--bai $INBAM.bai \
	-r $FASTA_FILE \
	--fai $FASTA_FILE.fai \
	--prefix $SAMPLE \
	--breakdancer --breakseq --manta --cnvnator --lumpy --genotype --svviz --delly_deletion --delly_insertion --delly_inversion --delly_duplication 

