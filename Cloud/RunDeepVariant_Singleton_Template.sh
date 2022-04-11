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
#SBATCH --array=0-29%6

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

# Pull latest version, if you already have it, this will be skipped
singularity pull docker://google/deepvariant:"${BIN_VERSION}"
export SINGULARITY_CACHEDIR=$PWD


# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

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

VCF=${SAMPLE}.vcf.gz
GVCF=${SAMPLE}.gvcf.gz

# Index CRAM
samtools index $BAM

# Run DeepVariant
singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
	-B "${BAM_DIR}":"/bamdir" \
	-B "${FASTA_DIR}":"/genomedir" \
	-B "${OUTPUT_DIR}":"/output" \
	-B "${TMP_DIR}":"/tmp/" \
	docker://google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --intermediate_results_dir="/tmp/" \
  --ref="/genomedir/$FASTA_FILE" \
  --reads="/bamdir/$SAMPLE.final.cram" \
  --output_vcf="/output/$VCF" \
  --output_gvcf="/output/$GVCF" \
  --num_shards=$NSLOTS 

