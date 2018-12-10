# New Pipeline
ID_FAMILY=<e.g. IM3>
ID_PROBAND=<e.g TIDEX-IM-008>
ID_MOTHER=<e.g. TIDEX-IM-009>
ID_FATHER=<e.g. TIDEX-IM-010>
PROBAND_SEX=<male or female>

DIR_SCRIPTS=/mnt/causes-vnx1/PIPELINES/AnnotateVariants/PipelineScripts/
SCRIPT_MASTER_WGS=${DIR_SCRIPTS}/Pipeline_Master_WGS.py

DIR_WORKING=/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/$ID_FAMILY/
DIR_RAW=/mnt/causes-vnx2/TIDE/RAW/EXOME_TIDEX/$ID_FAMILY/

READ_LENGTH=101


# Make a directory
mkdir -p $DIR_WORKING


# Initial Processing
python $SCRIPT_MASTER_WGS -v new \
 -d $DIR_WORKING -p 8 -m 30gb \
 -1 $DIR_RAW/${ID_PROBAND}_1.fastq.gz \
 -2 $DIR_RAW/${ID_PROBAND}_2.fastq.gz \
 -G GSC -s $ID_PROBAND --readlength $READ_LENGTH -S PBS

python $SCRIPT_MASTER_WGS -v new \
 -d $DIR_WORKING -p 8 -m 30gb \
 -1 $DIR_RAW/${ID_MOTHER}_1.fastq.gz \
 -2 $DIR_RAW/${ID_MOTHER}_2.fastq.gz \
 -G GSC -s $ID_MOTHER --readlength $READ_LENGTH -S PBS

python $SCRIPT_MASTER_WGS -v new \
 -d $DIR_WORKING -p 8 -m 30gb \
 -1 $DIR_RAW/${ID_FATHER}_1.fastq.gz \
 -2 $DIR_RAW/${ID_FATHER}_2.fastq.gz \
 -G GSC -s $ID_FATHER --readlength $READ_LENGTH -S PBS


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

# Change email address
# echo ""
# bash ~/SCRIPTS/replace_email_address.sh

