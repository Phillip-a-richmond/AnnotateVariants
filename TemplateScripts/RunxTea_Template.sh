#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
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
MINICONDA3_DIR=/mnt/common/Precision/Miniconda3/
source $MINICONDA3_DIR/opt/miniconda3/etc/profile.d/conda.sh
conda activate xTea


# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR
WORKING_DIR=/mnt/scratch/Precision/EPGEN/PROCESS/EPGEN083/

## Set working space
mkdir -p $WORKING_DIR
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=GRCh38
FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/
FASTA_FILE=GRCh38_full_analysis_set_plus_decoy_hla.fa
GENOME_FASTA=$FASTA_DIR/$FASTA_FILE

BAM_DIR=$WORKING_DIR
FAMILY_ID=EPGEN083
PROBAND_SAMPLEID=EPGEN083-01_${GENOME}
MOTHER_SAMPLEID=EPGEN083-02_${GENOME}
FATHER_SAMPLEID=EPGEN083-03_${GENOME}
SIBLING_SAMPLEID=sibling_id_${GENOME}
PED=$FAMILY_ID.ped

PROBAND_BAM=${PROBAND_SAMPLEID}.dupremoved.sorted.bam
FATHER_BAM=${FATHER_SAMPLEID}.dupremoved.sorted.bam
MOTHER_BAM=${MOTHER_SAMPLEID}.dupremoved.sorted.bam

# xTea stuff
ANALYSIS_DIR=$BAM_DIR/xTea
mkdir -p $ANALYSIS_DIR
XTEA_DIR=/mnt/common/Precision/xTea/xTea/xtea/
XTEA_ANNOTATION=/mnt/common/Precision/xTea/xTea/rep_lib_annotation/
XTEA_GFF=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GENCODE/gencode.v38.annotation.gff3



# ANNOTSV stuff
ANNOTSV=/mnt/common/Precision/AnnotSV/
export ANNOTSV=/mnt/common/Precision/AnnotSV/
GENELIST_BOOL=true
GENELIST=/mnt/scratch/Precision/EPGEN/PROCESS/EPGEN_Genes.txt

## Step 1 - Make Sample List
rm -f  ${FAMILY_ID}_sampleList.txt
echo "$PROBAND_SAMPLEID" >> ${FAMILY_ID}_sampleList.txt
echo "$FATHER_SAMPLEID" >> ${FAMILY_ID}_sampleList.txt
echo "$MOTHER_SAMPLEID" >> ${FAMILY_ID}_sampleList.txt

## Step 2 - Make file of listed alignments
rm -f ${FAMILY_ID}_sampleBamList.txt
echo "$PROBAND_SAMPLEID $BAM_DIR/$PROBAND_BAM" >> ${FAMILY_ID}_sampleBamList.txt
echo "$FATHER_SAMPLEID $BAM_DIR/$FATHER_BAM" >> ${FAMILY_ID}_sampleBamList.txt
echo "$MOTHER_SAMPLEID $BAM_DIR/$MOTHER_BAM" >> ${FAMILY_ID}_sampleBamList.txt

## Step 3 - Generate running script
xtea  \
	-i ${FAMILY_ID}_sampleList.txt \
	-b ${FAMILY_ID}_sampleBamList.txt \
	-x null \
	-p $ANALYSIS_DIR \
	-o xTea_Script_${FAMILY_ID}.sh \
	-l $XTEA_ANNOTATION \
	-r $GENOME_FASTA \
	-g $XTEA_GFF \
	--xtea $XTEA_DIR \
	-f 5907 \
	-y 15 \
	--slurm \
	-t 72:00:00 \
	-q dev_q \
	-m 40 \
	-n 5

# change email address:
for file in xTea/*/*/run_xTEA_pipeline.sh; do
	echo $file
	# add in a bit about the conda environment?
	sed -i 's/NONE/ALL/g' $file
	sed -i 's/chong.simonchu@gmail.com/prichmond@bcchr.ca/g' $file
done

exit


## Step 3 - Annotate
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment


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
