#!/bin/bash
  
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########

# Load environment 
ANNOTATE_VARIANTS_DIR=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATE_VARIANTS_DIR/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATE_VARIANTS_DIR/opt/AnnotateVariantsEnvironment/

# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Here working dir is the MELT subdirectory, assumes you have run MELT
WORKING_DIR=/mnt/scratch/Public/TESTING/GenomicsPipelineTest/MELT/

## Set working space
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=GRCh38
FASTA_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/
FASTA_FILE=GRCh38_full_analysis_set_plus_decoy_hla.fa
FAMILY_ID=EPGEN029

# BGZip the sub VCFs
bgzip ALU.final_comp.vcf
tabix ALU.final_comp.vcf.gz

bgzip HERVK.final_comp.vcf
tabix HERVK.final_comp.vcf.gz

bgzip SVA.final_comp.vcf
tabix SVA.final_comp.vcf.gz

bgzip LINE1.final_comp.vcf
tabix LINE1.final_comp.vcf.gz

# Concatenate the MELT files
bcftools concat -a \
       ALU.final_comp.vcf.gz \
       HERVK.final_comp.vcf.gz \
       SVA.final_comp.vcf.gz \
       LINE1.final_comp.vcf.gz \
	-o MergedMEI.vcf

export ANNOTSV=/mnt/common/Precision/AnnotSV/
# Annotate with targetd approach for epilepsy genes
$ANNOTSV/bin/AnnotSV -SVinputFile MergedMEI.vcf \
	-genomeBuild $GENOME \
        -candidateGenesFile /mnt/scratch/Precision/EPGEN/PROCESS/EPGEN_Genes.txt \
        -candidateGenesFiltering yes \
       	-outputFile ${FAMILY_ID}_MergedMEI.epilepsygenes.annotated.tsv 

$ANNOTSV/bin/AnnotSV -SVinputFile MergedMEI.vcf \
	-genomeBuild $GENOME \
       	-outputFile ${FAMILY_ID}_MergedMEI.annotated.tsv 



