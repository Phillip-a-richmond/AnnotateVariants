#!/bin/bash

#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=350G
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########

# Load env for bcftools/bedtools
ANNOTATEVARIANTS_INSTALL=annotate_variants_dir
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Load singularity
module load singularity

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

BAM_DIR=$WORKING_DIR
FAMILY_ID=family_id
PED=$FAMILY_ID.ped

PROBAND_BAM=${PROBAND_ID}_${GENOME}.dupremoved.sorted.bam
FATHER_BAM=${FATHER_ID}_${GENOME}.dupremoved.sorted.bam
MOTHER_BAM=${MOTHER_ID}_${GENOME}.dupremoved.sorted.bam

# Smooooooove 
# smoove - /mnt/common/Precision/Smoove/smoove_latest.sif
singularity run \
	-B /usr/lib/locale/:/usr/lib/locale/ \
	-B "${BAM_DIR}":"/bamdir" \
	-B "${FASTA_DIR}":"/genomedir" \
	-B "${OUTPUT_DIR}":"/output" \
	smoove_dir \
	smoove call \
	--outdir "/output/results-smoove/" \
	--duphold \
	-p $NSLOTS \
	-x --name ${FAMILY_ID} \
	--fasta "/genomedir/$FASTA_FILE" \
	--genotype \
	/bamdir/*dupremoved.sorted.bam \

# AnnotSV
export ANNOTSV=annotsv_dir
$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/results-smoove/${FAMILY_ID}-smoove.genotyped.vcf.gz \
	-genomeBuild genome_build \
	-overlap 80 \
	-reciprocal yes \
       	-outputFile ${FAMILY_ID}-smoove-annotsv 



