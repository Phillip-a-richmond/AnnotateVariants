# Pipeline Template December 20th, 2018 (The darkest day of the year!)

# Set variables
ID_FAMILY=STRAB
ID_PROBAND=STRAB-013
ID_MOTHER=STRAB-011
ID_FATHER=STRAB-012
PROBAND_SEX=female # male or female

DIR_SCRIPTS=/mnt/causes-vnx1/PIPELINES/AnnotateVariants/PipelineScripts/
SCRIPT_MASTER_WGS=${DIR_SCRIPTS}/Pipeline_Master_WGS.py
MTOOLBOX=/mnt/causes-vnx1/PIPELINES/AnnotateVariants/MToolBox_config_files/MToolBox_rCRS_config_with_markdup_and_indelrealign_RvdL.sh

DIR_WORKING=/mnt/causes-vnx2/TIDE/PROCESS/CYNTHIA_STRAB/
DIR_RAW=/mnt/causes-vnx2/TIDE/RAW/CYNTHIA_STRAB/

READ_LENGTH=101

EMAIL=prichmond@cmmt.ubc.ca


# Make a directory
mkdir -p $DIR_WORKING


# Initial Processing
python $SCRIPT_MASTER_WGS -v new \
	-d $DIR_WORKING -p 16 -m 30gb \
	-1 $DIR_RAW/${ID_PROBAND}_R1.fastq.gz \
	-2 $DIR_RAW/${ID_PROBAND}_R2.fastq.gz \
	-G GSC -s $ID_PROBAND --readlength $READ_LENGTH -S PBS --cnv \
	-E $EMAIL

python $SCRIPT_MASTER_WGS -v new \
	-d $DIR_WORKING -p 16 -m 30gb \
	-1 $DIR_RAW/${ID_MOTHER}_R1.fastq.gz \
	-2 $DIR_RAW/${ID_MOTHER}_R2.fastq.gz \
	-G GSC -s $ID_MOTHER --readlength $READ_LENGTH -S PBS --cnv \
	-E $EMAIL

python $SCRIPT_MASTER_WGS -v new \
	-d $DIR_WORKING -p 16 -m 30gb \
	-1 $DIR_RAW/${ID_FATHER}_R1.fastq.gz \
	-2 $DIR_RAW/${ID_FATHER}_R2.fastq.gz \
	-G GSC -s $ID_FATHER --readlength $READ_LENGTH -S PBS --cnv \
	-E $EMAIL


# Convert the gVCFs to GEMINI databases
python ${DIR_SCRIPTS}/Bam2Gemini.py -G GSC \
        -v GVCF -V ${ID_PROBAND}_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf,${ID_MOTHER}_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf,${ID_FATHER}_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf \
        -P $ID_FAMILY.ped -p 8 -m 30G -F $ID_FAMILY \
        -d $DIR_WORKING


# Make a PED file
python ${DIR_SCRIPTS}/MakePED.py \
        --proband ${ID_PROBAND}_BWAmem,${PROBAND_SEX},affected \
        --mother ${ID_MOTHER}_BWAmem,female,unaffected \
        --father ${ID_FATHER}_BWAmem,male,unaffected \
        --FamilyID $ID_FAMILY \
        -O $DIR_WORKING/$ID_FAMILY.ped

