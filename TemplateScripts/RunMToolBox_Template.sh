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

# Load environment 
ANNOTATE_VARIANTS_DIR=annotate_variants_dir
source $ANNOTATE_VARIANTS_DIR/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATE_VARIANTS_DIR/opt/AnnotateVariantsEnvironment/

# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR
WORKING_DIR=working_dir

## Set working space
mkdir -p $WORKING_DIR
cd $WORKING_DIR

FAMILY_ID=family_id
PED=$FAMILY_ID.ped
SAMPLE=sample_id
RAW_DIR=raw_dir
NSLOTS=$SLURM_CPUS_PER_TASK


FASTQR1=${RAW_DIR}/fastqr1
FASTQR2=${RAW_DIR}/fastqr2

echo $SAMPLE
ls $FASTQR1
ls $FASTQR2


# Run MToolBox for mitochondrial variant analysis 
echo "Mitochondrial Variant Analysis Started"
date
MTOOLBOX_PATH=mtoolbox_dir
PATH=$MTOOLBOX_PATH:$MTOOLBOX_PATH/MToolBox/:$PATH
MTOOLBOX_WORKING_DIR=$WORKING_DIR/MToolBox_${SAMPLE}/
mkdir -p $MTOOLBOX_WORKING_DIR/

MTOOLBOX_CONFIG_FILE_ORIGINAL='/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/MToolBox_config_files/MToolBox_RSRS_config_with_markdup_and_indelrealign_RvdL.sh'
MTOOLBOX_CONFIG_FILE_ORIGINAL_BASENAME=$(basename $MTOOLBOX_CONFIG_FILE_ORIGINAL)
MTOOLBOX_CONFIG_FILE=$MTOOLBOX_WORKING_DIR/$MTOOLBOX_CONFIG_FILE_ORIGINAL_BASENAME

#edit the MToolBox config file template so that is specifies the MToolBox results/working directory for the current analysis 
cp $MTOOLBOX_CONFIG_FILE_ORIGINAL $MTOOLBOX_CONFIG_FILE
sed -i "s#^output_name\=\.#output_name=$MTOOLBOX_WORKING_DIR#" $MTOOLBOX_CONFIG_FILE

#link the raw fastq files to the MToolBox working directory, and name them as required by MToolBox: \<sample\_name\>.R1.fastq, \<sample\_name\>.R2.fastq 
ln -sf ${FASTQR1} $MTOOLBOX_WORKING_DIR/${SAMPLE}.R1.fastq.gz
ln -sf ${FASTQR2} $MTOOLBOX_WORKING_DIR/${SAMPLE}.R2.fastq.gz

echo "Changing working directory to $MTOOLBOX_WORKING_DIR and running MToolBox..." 
PWD_CURRENT=`pwd`
cd $MTOOLBOX_WORKING_DIR
$MTOOLBOX_PATH/MToolBox/MToolBox.sh -m "-t $NSLOTS" -i ${MTOOLBOX_CONFIG_FILE}
echo "Changing working directory to back to $PWD_CURRENT..." 
cd $PWD_CURRENT



