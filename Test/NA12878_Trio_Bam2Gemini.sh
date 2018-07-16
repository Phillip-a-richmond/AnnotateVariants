#!/bin/bash
#PBS -N NA12878_Trio_Bam2Gemini
#PBS -V
#PBS -o /mnt/causes-vnx1/Pipelines/AnnotateVariants/Test/NA12878_Trio.o
#PBS -e /mnt/causes-vnx1/Pipelines/AnnotateVariants/Test/NA12878_Trio.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=40G
## Set the max walltime for the job
#PBS -l walltime=240:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=16
NSLOTS=$PBS_NUM_PPN
umask 0002
source /opt/tools/hpcenv.sh

FAMILY_ID='NA12878_Trio'
WORKING_DIR='/mnt/causes-vnx1/Pipelines/AnnotateVariants/Test/'
GENOME_FASTA='/mnt/causes-data01/data/GENOMES/hg19/FASTA/hg19.fa'
PED_FILE=$WORKING_DIR/NA12878_Trio.ped
TMPDIR=${WORKING_DIR}tmpdir/
mkdir $TMPDIR
SAMPLE1_VCF=NA12878_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf
SAMPLE2_VCF=NA12891_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf
SAMPLE3_VCF=NA12892_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf

# Step 2: Merge VCFs

#/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T CombineVariants \
#-R $GENOME_FASTA \
#--variant $WORKING_DIR${SAMPLE1_VCF} \
#--variant $WORKING_DIR${SAMPLE2_VCF} \
#--variant $WORKING_DIR${SAMPLE3_VCF} \
#-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf 
#
#Get Rid of non-chr chromosomes

#BGZIP that bad boy
#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf
#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz


# Step 3: Normalize merged VCF

# Define some variables

SNPEFFJAR=/opt/tools/snpEff/snpEff.jar
GEMINIDB=$WORKING_DIR${FAMILY_ID}.db
VCF=$WORKING_DIR${FAMILY_ID}.merged.hc.vcf
NORMVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcf.gz
zless $VCF \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| /opt/tools/vt/vt decompose -s - \
	| /opt/tools/vt/vt normalize -r $GENOME_FASTA - \
	| java -Xmx10g -jar $SNPEFFJAR GRCh37.75 \
	| /opt/tools/tabix/bgzip -c > $NORMVCF 
/opt/tools/tabix/tabix -p vcf $NORMVCF

## Step 4: Filter Merged, normalized VCF

NORMFILTERVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.filter.vcf.gz
/opt/tools/bcftools-1.8/bin/bcftools filter \
	 --include 'FORMAT/AD[0:1]>=10 && FORMAT/DP[0] < 300' \
	 -m + \
	 -s + \
	 -O z \
	 --output $NORMFILTERVCF \
	 $NORMVCF 

/opt/tools/tabix/tabix $NORMFILTERVCF \


# Step 5: VCFAnno - Turn your VCF file into an annotated VCF file
ANNOVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcfanno.vcf.gz 
/opt/tools/vcfanno/vcfanno -lua /mnt/causes-data01/data/RICHMOND/AnnotateVariants/VCFAnno/custom.lua \
-p $NSLOTS \
/mnt/causes-vnx1/Pipelines/AnnotateVariants/VCFAnno/VCFANNO_Config_PlusGNOMAD_PlusInHouse_SplitByPop_gnomAD_Exome.toml \
$NORMFILTERVCF > $ANNOVCF 


# Step 6: VCF2DB - Turn your annotated VCF file into a GEMINI DB

python /opt/tools/vcf2db/vcf2db.py \
--expand gt_quals --expand gt_depths --expand gt_alt_depths --expand gt_ref_depths --expand gt_types \
 --a-ok InHouseDB_AC  --a-ok in_segdup --a-ok AF --a-ok AC --a-ok AN --a-ok MLEAC --a-ok MLEAF --a-ok gnomad_genome_hom_all --a-ok gnomad_genome_hom_afr --a-ok gnomad_genome_hom_amr --a-ok gnomad_genome_hom_asj --a-ok gnomad_genome_hom_eas --a-ok gnomad_genome_hom_fin --a-ok gnomad_genome_hom_nfe --a-ok gnomad_genome_hom_oth --a-ok gnomad_exome_hom_all --a-ok gnomad_exome_hom_afr --a-ok gnomad_exome_hom_amr --a-ok gnomad_exome_hom_asj --a-ok gnomad_exome_hom_eas --a-ok gnomad_exome_hom_fin --a-ok gnomad_exome_hom_nfe --a-ok gnomad_exome_hom_oth --a-ok num_exac_Het --a-ok num_exac_Hom --a-ok cpg_island --a-ok common_pathogenic --a-ok cse-hiseq --a-ok DS --a-ok ConfidentRegion \
$ANNOVCF $PED_FILE $GEMINIDB 
