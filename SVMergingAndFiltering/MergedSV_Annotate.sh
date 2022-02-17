#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
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

#################
# For Merged SV #
#################


##########
# Smoove #
##########
#INPUT_MERGED_VCF=InHouseDB_SV_50_20211130.vcf
#INPUT_MERGED_VCF_CLEANED=InHouseDB_SV_50_20211130_noTRA.vcf
#INPUT_MERGED_VCF_BND=InHouseDB_SV_50_20211130_BND.vcf
#OUTPUT=InHouseDB_SV_50_20211130_noTRA_annotsv
#
## Clean out TRA (from BNDs)
#grep -v -w 'TRA' $INPUT_MERGED_VCF > $INPUT_MERGED_VCF_CLEANED
#
#
## AnnotSV
### If GeneList
#if [ "$GENELIST_BOOL" = true ]; then
#	$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/$INPUT_MERGED_VCF_CLEANED \
#		-genomeBuild $GENOME \
#	        -candidateGenesFile $GENELIST \
#                -candidateGenesFiltering yes \
#	       	-outputFile ${OUTPUT}-candidateGenes
#fi
#
#
## Normal
#$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/$INPUT_MERGED_VCF_CLEANED \
#	-genomeBuild $GENOME \
#       	-outputFile $OUTPUT


# Focus on BNDs
# In development
#grep 'SVTYPE=TRA' $INPUT_MERGED_VCF | sed -e 's/TRA/BND/g' > InHouseDB_SV_50_20211130_BND.vcf



#########
# MANTA #
#########
INPUT_MERGED_VCF=InHouseDB_SV_Manta_50_20211130.vcf
INPUT_MERGED_VCF_CLEANED=InHouseDB_SV_Manta_50_20211130_noTRA.vcf
OUTPUT=InHouseDB_SV_Manta_50_20211130_noTRA_annotsv

# Clean out TRA (from BNDs)
grep -v -w 'TRA' $INPUT_MERGED_VCF > $INPUT_MERGED_VCF_CLEANED


# AnnotSV
## If GeneList
if [ "$GENELIST_BOOL" = true ]; then
	$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/$INPUT_MERGED_VCF_CLEANED \
		-genomeBuild $GENOME \
	        -candidateGenesFile $GENELIST \
                -candidateGenesFiltering yes \
	       	-outputFile ${OUTPUT}-candidateGenes
fi


# Normal
$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/$INPUT_MERGED_VCF_CLEANED \
	-genomeBuild $GENOME \
       	-outputFile $OUTPUT

exit

# Focus on BNDs
# In development
#grep 'SVTYPE=TRA' $INPUT_MERGED_VCF | sed -e 's/TRA/BND/g' > InHouseDB_SV_50_20211130_BND.vcf


##################
# For Merged MEI #
##################
INPUT_MERGED_VCF=InHouseDB_MEI_15_20211130.vcf
OUTPUT=InHouseDB_MEI_15_20211130_annotsv
## If GeneList
if [ "$GENELIST_BOOL" = true ]; then
        $ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/$INPUT_MERGED_VCF \
                -genomeBuild $GENOME \
                -candidateGenesFile $GENELIST \
                -candidateGenesFiltering yes \
                -outputFile ${OUTPUT}-candidateGenes 
fi


# Normal
$ANNOTSV/bin/AnnotSV -SVinputFile $WORKING_DIR/$INPUT_MERGED_VCF \
        -genomeBuild $GENOME \
        -outputFile $OUTPUT


