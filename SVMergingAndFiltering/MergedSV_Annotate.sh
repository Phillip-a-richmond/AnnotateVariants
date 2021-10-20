#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
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
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Smoove
# Assumes you have a smoove.sh script copied to your BAM dir
SMOOVE_SIF=/mnt/common/Precision/Smoove/smoove_latest.sif
ANNOTSV=/mnt/common/Precision/AnnotSV/

# Load singularity
module load singularity

# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR
WORKING_DIR=/mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/

## Set working space
mkdir -p $WORKING_DIR
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=GRCh38
FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/
FASTA_FILE=GRCh38_full_analysis_set_plus_decoy_hla.fa

BAM_DIR=$WORKING_DIR

# ANNOTSV Stuff
GENELIST_BOOL=true
GENELIST=/mnt/scratch/Precision/EPGEN/PROCESS/EPGEN_Genes.txt

export ANNOTSV=/mnt/common/Precision/AnnotSV/
# AnnotSV
# For Merged SV
## If GeneList
if [ "$GENELIST_BOOL" = true ]; then
	$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/InHouseDB_SV_50_20210827_noTRA.vcf \
		-genomeBuild $GENOME \
	        -candidateGenesFile $GENELIST \
                -candidateGenesFiltering yes \
	       	-outputFile InHouseDB_SV_50_20210827_noTRA-annotsv-candidateGenes 
fi


## Normal
#$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/InHouseDB_SV_50_20210827_noTRA.vcf \
#	-genomeBuild $GENOME \
#       	-outputFile InHouseDB_SV_50_20210827_noTRA-annotsv 


# For Merged MEI
# AnnotSV
# For Merged SV
## If GeneList
if [ "$GENELIST_BOOL" = true ]; then
        $ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/InHouseDB_MEI_15_20210827.vcf \
                -genomeBuild $GENOME \
                -candidateGenesFile $GENELIST \
                -candidateGenesFiltering yes \
                -outputFile InHouseDB_MEI_15_20210827-annotsv-candidateGenes
fi


## Normal
#$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/InHouseDB_MEI_15_20210827.vcf \
#        -genomeBuild $GENOME \
#        -outputFile InHouseDB_MEI_15_20210827-annotsv


