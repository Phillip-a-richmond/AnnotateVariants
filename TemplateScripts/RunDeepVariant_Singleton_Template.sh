#!/bin/bash

#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=350G
#SBATCH --cpus-per-task=40
#SBATCH --time=72:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########

# Load singularity
module load singularity

# Point to the deepvariant executable
DeepVariant_SIF=deepvariant_sif

# Point to the GLNexus_CLI
GLNexus_CLI=glnexus_cli

# Load env for bcftools
ANNOTATEVARIANTS_INSTALL=annotate_variants_dir
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Singularity cache
unset $PYTHONPATH
export SINGULARITY_CACHEDIR=$PWD

# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR
WORKING_DIR=working_dir
OUTPUT_DIR=$WORKING_DIR

## Set working space
mkdir -p $WORKING_DIR
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=genome_build
FASTA_DIR=fasta_dir
FASTA_FILE=fasta_file

SEQ_TYPE=seq_type
BAM_DIR=$WORKING_DIR
FAMILY_ID=family_id
PROBAND_ID=proband_id
PED=$FAMILY_ID.ped


PROBAND_BAM=${PROBAND_ID}_${GENOME}.dupremoved.sorted.bam

PROBAND_VCF=${PROBAND_ID}.vcf.gz

PROBAND_GVCF=${PROBAND_ID}.gvcf.gz

# Proband 
singularity exec -e -c -B /usr/lib/locale/:/usr/lib/locale/ \
	-B "${BAM_DIR}":"/bamdir" \
	-B "${FASTA_DIR}":"/genomedir" \
	-B "${OUTPUT_DIR}":"/output" \
	-W $OUTPUT_DIR \
	$DeepVariant_SIF \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/genomedir/$FASTA_FILE" \
  --intermediate_results_dir="/output/intermediate_results_dir" \
  --reads="/bamdir/$PROBAND_BAM" \
  --output_vcf="/output/$PROBAND_VCF" \
  --output_gvcf="/output/$PROBAND_GVCF" \
  --num_shards=$NSLOTS 

cp $PROBAND_VCF ${FAMILY_ID}.merged.vcf.gz

