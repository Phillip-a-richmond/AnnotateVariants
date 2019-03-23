#!/bin/bash
#PBS -N NA12878_Trio_Bam2Gemini
#PBS -V
#PBS -o /mnt/causes-vnx1/PIPELINES/AnnotateVariants/Test/NA12878_Trio.o
#PBS -e /mnt/causes-vnx1/PIPELINES/AnnotateVariants/Test/NA12878_Trio.e
## Set the total memory for the job
#PBS -l mem=40G
## Set the max walltime for the job
#PBS -l walltime=240:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=16
NSLOTS=1
umask 0002
source /opt/tools/hpcenv.sh

# Define Tool paths. If they are in your path, simply change these full filepaths to only be the final command
# For example: Change BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools to be BCFTOOLS=bcftools if it's in your path 

ANNOTVARDIR=/mnt/causes-vnx1/PIPELINES/AnnotateVariants/
SNPEFFJAR=/opt/tools/snpEff/snpEff.jar
BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools
VCFANNO=/opt/tools/vcfanno/vcfanno
VCF2DB=/opt/tools/vcf2db/vcf2db.py
GATKJAR=/opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar
JAVA=/opt/tools/jdk1.7.0_79/bin/java
BGZIP=/opt/tools/tabix/bgzip
TABIX=/opt/tools/tabix/tabix
VT=/opt/tools/vt/vt

# Define variables for what you are working with

FAMILY_ID='NA12878_Trio'
WORKING_DIR='/mnt/causes-vnx1/PIPELINES/AnnotateVariants/Test/'
GENOME_FASTA='/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa'
PED_FILE=$WORKING_DIR/NA12878_Trio.ped
TMPDIR=${WORKING_DIR}tmpdir/
mkdir $TMPDIR
GEMINIDB=$WORKING_DIR${FAMILY_ID}.db
VCF=$WORKING_DIR${FAMILY_ID}.merged.hc.vcf
NORMVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcf.gz
NORMFILTERVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.filter.vcf.gz
ANNOVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcfanno.vcf.gz 
SAMPLE1_VCF=NA12878_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf
SAMPLE2_VCF=NA12891_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf
SAMPLE3_VCF=NA12892_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf

# Merge VCFs

#$JAVA -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -T CombineVariants \
#-R $GENOME_FASTA \
#--variant $WORKING_DIR${SAMPLE1_VCF} \
#--variant $WORKING_DIR${SAMPLE2_VCF} \
#--variant $WORKING_DIR${SAMPLE3_VCF} \
#-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf 
#
##Get Rid of non-chr chromosomes
#
##  Normalize merged VCF, annotate with SNPeff
#
#zless $VCF \
#	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
#	| $VT decompose -s - \
#	| $VT normalize -r $GENOME_FASTA - \
#	| $JAVA -Xmx10g -jar $SNPEFFJAR GRCh37.75 \
#	| $BGZIP -c > $NORMVCF 
#$TABIX -p vcf $NORMVCF
#
## Filter Merged, normalized VCF
#
#$BCFTOOLS filter \
#	 --include 'FORMAT/AD[*:1]>=5 && FORMAT/DP[*] < 600' \
#	 -m + \
#	 -s + \
#	 -O z \
#	 --output $NORMFILTERVCF \
#	 $NORMVCF 
#
#$TABIX $NORMFILTERVCF \


# VCFAnno - Turn your VCF file into an annotated VCF file
$VCFANNO -lua $ANNOTVARDIR/VCFAnno/custom.lua \
-p 16 -base-path /mnt/causes-vnx1/DATABASES/ \
$ANNOTVARDIR/VCFAnno/VCFAnno_Config_20190321_GAC.toml \
$NORMFILTERVCF > $ANNOVCF 


# VCF2DB - Turn your annotated VCF file into a GEMINI DB

python $VCF2DB \
 --expand gt_quals --expand gt_depths --expand gt_alt_depths --expand gt_ref_depths --expand gt_types \
 --a-ok gnomad_exome_ac_global --a-ok gnomad_exome_ac_popmax --a-ok gnomad_exome_an_global --a-ok gnomad_exome_an_popmax --a-ok gnomad_exome_hom_controls --a-ok gnomad_exome_hom_global \
 --a-ok gnomad_exome_hom_popmax --a-ok gnomad_exome_popmax --a-ok gnomad_genome_ac_global --a-ok gnomad_genome_ac_popmax --a-ok gnomad_genome_an_global --a-ok gnomad_genome_an_popmax \
 --a-ok gnomad_genome_hom_controls --a-ok gnomad_genome_hom_global --a-ok gnomad_genome_hom_popmax --a-ok gnomad_genome_popmax \
 --a-ok InHouseDB_AC  --a-ok in_segdup --a-ok AF --a-ok AC --a-ok AN --a-ok MLEAC --a-ok MLEAF --a-ok cpg_island --a-ok common_pathogenic --a-ok cse-hiseq --a-ok DS --a-ok ConfidentRegion \
$ANNOVCF $PED_FILE $GEMINIDB 
