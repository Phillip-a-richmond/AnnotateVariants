#!/bin/bash
  
#SBATCH --mail-user=email_address
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

# Load env for bcftools
ANNOTATEVARIANTS_INSTALL=annotate_variants_dir
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment
export ANNOTSV=annotsv_dir

# Here working dir is the MELT subdirectory, assumes you have run MELT
WORKING_DIR=working_dir/MELT/

## Set working space
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=genome_build
FASTA_DIR=fasta_dir
FASTA_FILE=fasta_file

# ANNOTSV Stuff
GENELIST_BOOL=genelist_bool
GENELIST=genelist

FAMILY_ID=family_id

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

# Annotate with targetd approach for epilepsy genes
if [ "$GENELIST_BOOL" = true ]; then
	$ANNOTSV/bin/AnnotSV -SVinputFile MergedMEI.vcf \
		-genomeBuild $GENOME \
	        -candidateGenesFile $GENELIST \
	        -candidateGenesFiltering yes \
	       	-outputFile ${FAMILY_ID}_MergedMEI-candidateGenes-annotsv 
fi

$ANNOTSV/bin/AnnotSV -SVinputFile MergedMEI.vcf \
	-genomeBuild $GENOME \
       	-outputFile ${FAMILY_ID}_MergedMEI-annotsv



