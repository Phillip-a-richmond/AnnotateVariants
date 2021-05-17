#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
## 160 Gb of RAM for the whole job
#SBATCH --mem=160G

## Using 40 CPUs
#SBATCH --cpus-per-task=40

## Running for a max time of 48 hours
#SBATCH --time=48:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

NSLOTS=$SLURM_CPUS_PER_TASK

# Load the necessary tools
source /mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/opt/miniconda3/etc/profile.d/conda.sh 
conda activate /mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/opt/AnnotateVariantsEnvironment/

# Software
ANNOTVARDIR=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
BCFTOOLS=bcftools
VCFANNO=vcfanno
VCF2DB=vcf2db.py
JAVA=java
BGZIP=bgzip
TABIX=tabix
# Databases
GENOME_FASTA=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa
LOFSCORESFILE=/mnt/common/DATABASES/GENERIC/VEP/LoFtool_scores.txt
SPLICEAISNV=/mnt/common/DATABASES/REFERENCES/GRCh38/SPLICEAI/spliceai_scores.masked.snv.hg38.vcf.gz
SPLICEAIINDEL=/mnt/common/DATABASES/REFERENCES/GRCh38/SPLICEAI/spliceai_scores.masked.indel.hg38.vcf.gz
CACHEDIR=/mnt/common/DATABASES/REFERENCES/GRCh38/VEP/
PLUGINSDIR=/mnt/common/DATABASES/REFERENCES/GRCh38/VEP/PLUGINS/
MAXENTSCANDIR=/mnt/common/WASSERMAN_SOFTWARE/VEP/fordownload/
SLIVAR=/mnt/common/Precision/Slivar/
gff=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/Homo_sapiens.GRCh38.100.chr.gff3.gz

# Files
WORKDIR=/mnt/scratch/Precision/EPGEN/PROCESS/EPGEN012_PAR/
cd $WORKDIR

SAMPLE=EPGEN012
PED=EPGEN012.ped
INVCF=${SAMPLE}.glnexus.merged.vcf.gz

# check for VCF
if [ ! -f $INVCF ]; then
	echo "Missing vcf"
	exit
fi


# If you find it, echo the name
echo "VCF found"
ls -lahtr $INVCF

INVCF_NORM=${SAMPLE}_bcftoolsnorm.vcf.gz
INVCF_NORM_FILTER=${SAMPLE}_bcftoolsnorm.bcftoolsfilter.vcf.gz
INVCF_NORM_FILTER_NOREFCALL=${SAMPLE}_bcftoolsnorm.bcftoolsfilter.NoRefCall.vcf.gz
INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE=${SAMPLE}_bcftoolsnorm.bcftoolsfilter.NoRefCall.slivarrare.vcf.gz
INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO=${SAMPLE}_bcftoolsnorm.bcftoolsfilter.NoRefCall.slivarrare.vepanno.vcf.gz
INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT=${SAMPLE}_bcftoolsnorm.bcftoolsfilter.NoRefCall.slivarrare.vepanno.vepsplit.vcf.gz
INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT_VCFANNO=${SAMPLE}_bcftoolsnorm.bcftoolsfilter.NoRefCall.slivarrare.vepanno.vepsplit.vcfanno.vcf
INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT_VCFANNO_NUMBERFIX=${SAMPLE}_bcftoolsnorm.bcftoolsfilter.NoRefCall.slivarrare.vepanno.vepsplit.vcfanno.numberfix.vcf
GEMINIDB=${SAMPLE}_Gemini.db

# Cutoffs
HOM_ALT_CUTOFF=15
AF_CUTOFF=0.001
# Step 1 - BCFTools normalize & apply filter
# BCFTOOLS
# Normalize VCF
$BCFTOOLS norm \
	-f $GENOME_FASTA \
	--threads @NSLOTS \
	-m - \
	-O z \
	--output $INVCF_NORM \
	$INVCF 

$TABIX -f $INVCF_NORM

# Filter normalized VCF for high depth + low alt allele support (soft-filter)
$BCFTOOLS filter \
         --include 'FORMAT/AD[*:1]>=5 && FORMAT/DP[*] < 600' \
         -m + \
         -s + \
         -O z \
         --output $INVCF_NORM_FILTER \
         $INVCF_NORM

$TABIX -f $INVCF_NORM_FILTER

# Step 2 - Get rid of RefCall from DeepVariant output
zgrep -v "RefCall" $INVCF_NORM_FILTER | bgzip -c > $INVCF_NORM_FILTER_NOREFCALL
tabix $INVCF_NORM_FILTER_NOREFCALL


# Step 3 - Slivar filter for gnomad frequency, here set to GRCh38
#SLIVAR 
$SLIVAR/slivar.exe expr --vcf $INVCF_NORM_FILTER_NOREFCALL \
    --ped $PED \
    --pass-only \
    -g $SLIVAR/gnomad.hg38.genomes.v3.fix.zip \
    --info "INFO.gnomad_popmax_af < $AF_CUTOFF && INFO.gnomad_nhomalt < $HOM_ALT_CUTOFF && variant.FILTER == \"PASS\" && variant.ALT[0] != \"*\"" \
    --js $SLIVAR/slivar/js/slivar-functions.js \
    -o $INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE

# Step 4 - Run VEP to annotate
# VEP
vep \
	-i $INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE \
	--buffer_size 10000 \
	--polyphen b\
	--fasta $GENOME_FASTA \
	--sift b \
	--biotype \
	--hgvs \
	--protein \
	--domains \
	--pubmed \
	--ccds \
	--check_existing \
	--nearest symbol \
	--gene_phenotype \
	--canonical \
	--force_overwrite \
	--offline \
	--cache \
	--dir_cache $CACHEDIR \
	--dir_plugins $PLUGINSDIR \
	--plugin MaxEntScan,$MAXENTSCANDIR \
	--plugin SpliceAI,snv=$SPLICEAISNV,indel=$SPLICEAIINDEL \
	-o $INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO \
	--compress_output bgzip \
	--fork 20 \
	--vcf 


# Step 5 - Split VEP annotations into unique records
# BCFTOOLS
bcftools +split-vep \
	-p vep \
	-a CSQ \
	-O z \
	-o $INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT \
	-c MaxEntScan_alt:Float,\
MaxEntScan_diff:Float,\
MaxEntScan_ref:Float,\
SpliceAI_pred_DP_AG:Float,\
SpliceAI_pred_DP_AL:Float,\
SpliceAI_pred_DP_DG:Float,\
SpliceAI_pred_DP_DL:Float,\
SpliceAI_pred_DS_AG:Float,\
SpliceAI_pred_DS_AL:Float,\
SpliceAI_pred_DS_DG:Float,\
SpliceAI_pred_DS_DL:Float,\
SpliceAI_pred_SYMBOL:String,\
Existing_variation:String,\
Protein_position:String,\
DOMAINS:String,\
PUBMED:String,\
CLIN_SIG:String,\
	-s worst \
	$INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO


# Step 6 - VCFAnno - Turn your VCF file into an annotated VCF file
# VCFANNO
$VCFANNO -lua $ANNOTVARDIR/VCFAnno/custom.lua \
-p $NSLOTS -base-path /mnt/common/DATABASES/REFERENCES/GRCh38/ \
$ANNOTVARDIR/VCFAnno/VCFAnno_Config_20210222_GPCC_NoSpliceAI_GRCh38.toml \
$INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT > $INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT_VCFANNO

# Step 7 - Fix header for VEP-added annotations
## Fix the Number=1 --> Number=.
## This is necessary to read these as floats:
##INFO=<ID=vepMaxEntScan_alt,Number=.,Type=Float,Description="The MaxEntScan_alt field from INFO/CSQ">
##INFO=<ID=vepMaxEntScan_diff,Number=.,Type=Float,Description="The MaxEntScan_diff field from INFO/CSQ">
##INFO=<ID=vepMaxEntScan_ref,Number=.,Type=Float,Description="The MaxEntScan_ref field from INFO/CSQ">
##INFO=<ID=vepSpliceAI_pred_DP_AG,Number=.,Type=Float,Description="The SpliceAI_pred_DP_AG field from INFO/CSQ">
##INFO=<ID=vepSpliceAI_pred_DP_AL,Number=.,Type=Float,Description="The SpliceAI_pred_DP_AL field from INFO/CSQ">
##INFO=<ID=vepSpliceAI_pred_DP_DG,Number=.,Type=Float,Description="The SpliceAI_pred_DP_DG field from INFO/CSQ">
##INFO=<ID=vepSpliceAI_pred_DP_DL,Number=.,Type=Float,Description="The SpliceAI_pred_DP_DL field from INFO/CSQ">
##INFO=<ID=vepSpliceAI_pred_DS_AG,Number=.,Type=Float,Description="The SpliceAI_pred_DS_AG field from INFO/CSQ">
##INFO=<ID=vepSpliceAI_pred_DS_AL,Number=.,Type=Float,Description="The SpliceAI_pred_DS_AL field from INFO/CSQ">
##INFO=<ID=vepSpliceAI_pred_DS_DG,Number=.,Type=Float,Description="The SpliceAI_pred_DS_DG field from INFO/CSQ">
##INFO=<ID=vepSpliceAI_pred_DS_DL,Number=.,Type=Float,Description="The SpliceAI_pred_DS_DL field from INFO/CSQ">

sed -e 's/vepMaxEntScan_alt,Number=./vepMaxEntScan_alt,Number=1/g' \
	-e 's/vepMaxEntScan_diff,Number=./vepMaxEntScan_diff,Number=1/g'  \
	 -e 's/vepMaxEntScan_ref,Number=./vepMaxEntScan_ref,Number=1/g'   \
	 -e 's/vepSpliceAI_pred_DP_AG,Number=./vepSpliceAI_pred_DP_AG,Number=1/g'  \
	 -e 's/vepSpliceAI_pred_DP_AL,Number=./vepSpliceAI_pred_DP_AL,Number=1/g'  \
	 -e 's/vepSpliceAI_pred_DP_DG,Number=./vepSpliceAI_pred_DP_DG,Number=1/g'  \
	 -e 's/vepSpliceAI_pred_DP_DL,Number=./vepSpliceAI_pred_DP_DL,Number=1/g'  \
	 -e 's/vepSpliceAI_pred_DS_AG,Number=./vepSpliceAI_pred_DS_AG,Number=1/g'  \
	 -e 's/vepSpliceAI_pred_DS_AL,Number=./vepSpliceAI_pred_DS_AL,Number=1/g'  \
	 -e 's/vepSpliceAI_pred_DS_DG,Number=./vepSpliceAI_pred_DS_DG,Number=1/g'  \
	 -e 's/vepSpliceAI_pred_DS_DL,Number=./vepSpliceAI_pred_DS_DL,Number=1/g'  \
	$INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT_VCFANNO  \
	> $INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT_VCFANNO_NUMBERFIX


# Step 8 - Create a GEMINI-compatible database
# VCF2DB

### Cover these for --a-ok (anything with number=A: 
### ##INFO=<ID=gnomad_genome_ac_global,Number=A,Type=Integer,Description="Alternate allele count for samples (from /mnt/common/DATABASES/REFERENCES/GRCh38//GNOMAD/V3/gnomad.genomes.r3.0.sites.vcf.gz)">
####INFO=<ID=gnomad_genome_af_global,Number=A,Type=Float,Description="Alternate allele frequency in samples (from /mnt/common/DATABASES/REFERENCES/GRCh38//GNOMAD/V3/gnomad.genomes.r3.0.sites.vcf.gz)">
###INFO=<ID=gnomad_genome_hom_global,Number=A,Type=Integer,Description="Count of homozygous individuals in samples (from /mnt/common/DATABASES/REFERENCES/GRCh38//GNOMAD/V3/gnomad.genomes.r3.0.sites.vcf.gz)">
###FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">

$VCF2DB \
 --expand gt_quals --expand gt_depths --expand gt_alt_depths --expand gt_ref_depths --expand gt_types \
 --a-ok AF --a-ok AC --a-ok AN \
 --a-ok gnomad_genome_ac_global --a-ok gnomad_genome_af_global --a-ok gnomad_genome_hom_global --a-ok VAF \
 --a-ok cosmic_count_observed --a-ok gnomad_genome_an_global --a-ok gnomad_nhomalt \
$INVCF_NORM_FILTER_NOREFCALL_SLIVARRARE_VEPANNO_SPLIT_VCFANNO_NUMBERFIX $PED $GEMINIDB 


