#!/bin/bash

#SBATCH --account=rrg-wyeth

## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --time=2-0:00
#SBATCH --nodes=1
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


WORKING_DIR=/scratch/richmonp/DATABASES/CLINVAR/
VCF=Clinvar_norm_eff.vcf
ANNOVCF=Clinvar_norm_eff_anno.vcf
TMPDIR=$WORKING_DIR/tmp
ANNOTATE_VARIANTS_DIR=/home/richmonp/scratch/AnnotateVariants/
rm -rf $TMPDIR
mkdir $TMPDIR

cd $WORKING_DIR

vcfanno -p 8 \
	-lua ${ANNOTATE_VARIANTS_DIR}rare-disease.lua \
	${ANNOTATE_VARIANTS_DIR}VCFANNO_Config_Cedar.toml \
	$VCF > $ANNOVCF.3 


