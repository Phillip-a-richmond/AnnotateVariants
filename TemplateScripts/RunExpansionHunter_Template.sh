#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL
#SBATCH -p dev_q

## CPU Usage
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########

# Run ExpansionHunter DeNovo 
EH5=e5h_var
# eh5/mnt/common/Precision/ExpansionHunter/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter
CATALOG=eh5_catalog_var
# /mnt/common/Precision/ExpansionHunter/ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog/hg38/variant_catalog.json

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
SIBLING_BAM=${SIBLING_SAMPLEID}.dupremoved.sorted.bam

PROBAND_BASE=${PROBAND_SAMPLEID}_EH
FATHER_BASE=${FATHER_SAMPLEID}_EH
MOTHER_BASE=${MOTHER_SAMPLEID}_EH
SIBLING_BASE=${SIBLING_SAMPLEID}_EH

MOTHER_PRESENT=mother_boolean
FATHER_PRESENT=father_boolean
SIBLING_PRESENT=sibling_boolean


## Step 1 - Run EH on BAM files
### Proband
$EH5 --reference $FASTA_DIR/$FASTA_FILE \
	--variant-catalog $CATALOG \
	--output-prefix $PROBAND_BASE \
	--threads $NSLOTS \
	--reads $PROBAND_BAM

### Father
if [ "$FATHER_PRESENT" = true ] ; then
$EH5 --reference $FASTA_DIR/$FASTA_FILE \
        --variant-catalog $CATALOG \
        --output-prefix $FATHER_BASE \
        --threads $NSLOTS \
        --reads $FATHER_BAM
fi

### Mother
if [ "$MOTHER_PRESENT" = true ] ; then
$EH5 --reference $FASTA_DIR/$FASTA_FILE \
        --variant-catalog $CATALOG \
        --output-prefix $MOTHER_BASE \
        --threads $NSLOTS \
        --reads $MOTHER_BAM
fi

### Sibling
if [ "$SIBLING_PRESENT" = true ] ; then
$EH5 --reference $FASTA_DIR/$FASTA_FILE \
	--variant-catalog $CATALOG \
	--output-prefix $SIBLING_BASE \
	--threads $NSLOTS \
	--reads $SIBLING_BAM
fi


