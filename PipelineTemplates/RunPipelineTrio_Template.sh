#############
### SETUP ###
#############
## FAMILY INFORMATION AND SAMPLE IDs
ID_FAMILY=<e.g. IM3>
ID_PROBAND=<e.g TIDEX-IM-008>
ID_MOTHER=<e.g. TIDEX-IM-009>
ID_FATHER=<e.g. TIDEX-IM-010>

PROBAND_SEX=<male or female>


## EMAIL ADDRESS FOR JOB NOTIFICATIONS
EMAIL=<e.g. rvdlee@cmmt.ubc.ca>


## DIRECTORIES AND SCRIPT PATHS
BASE_DIR=/mnt/causes-vnx1/PIPELINES/AnnotateVariants/
DIR_SCRIPTS=${BASE_DIR}/PipelineScripts/
DB_DIR=/mnt/causes-vnx1/DATABASES/

SCRIPT_MASTER=${DIR_SCRIPTS}/Pipeline_Master.py
SCRIPT_BAM2GEMINI=${DIR_SCRIPTS}/Bam2Gemini.py
PEDMAKER=${BASE_DIR}PipelineScripts/MakePED.py

MTOOLBOX_CONFIG_FILE=${BASE_DIR}/MToolBox_config_files/MToolBox_rCRS_config_with_markdup_and_indelrealign_RvdL.sh

DIR_RAW=/mnt/causes-vnx2/TIDE/RAW/EXOME_TIDEX/$ID_FAMILY/
DIR_WORKING=/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/$ID_FAMILY/

READ_LENGTH=151 # Obtain information by running:    $DIR_SCRIPTS/read_length.sh $DIR_RAW/${ID_PROBAND}_1.fastq.gz (and other IDs, _1, _2 etc). Or for all *fast.gz:    find $DIR_RAW/ -name "*.gz" -exec $DIR_SCRIPTS/read_length.sh {} \;

mkdir -p $DIR_WORKING # Make a directory


###############################
### CREATE PIPELINE SCRIPTS ###
###############################
## STEP 1) processing fastq to bam to vcf, creates "FullPipeline" scripts

# proband
python $SCRIPT_MASTER -v new \
		 -T Exome \
		 -d $DIR_WORKING -p 12 -m 40gb \
		 -1 $DIR_RAW/${ID_PROBAND}_1.fastq.gz \
		 -2 $DIR_RAW/${ID_PROBAND}_2.fastq.gz \
		 -E $EMAIL \
		 -G GSC -s $ID_PROBAND --readlength $READ_LENGTH -S PBS \
		 --mtoolbox $MTOOLBOX_CONFIG_FILE --metrics-exome

# mother
python $SCRIPT_MASTER -v new \
		 -T Exome \
		 -d $DIR_WORKING -p 12 -m 40gb \
		 -1 $DIR_RAW/${ID_MOTHER}_1.fastq.gz \
		 -2 $DIR_RAW/${ID_MOTHER}_2.fastq.gz \
		 -E $EMAIL \
		 -G GSC -s $ID_MOTHER --readlength $READ_LENGTH -S PBS \
		 --mtoolbox $MTOOLBOX_CONFIG_FILE --metrics-exome

# father
python $SCRIPT_MASTER -v new \
		 -T Exome \
		 -d $DIR_WORKING -p 12 -m 40gb \
		 -1 $DIR_RAW/${ID_FATHER}_1.fastq.gz \
		 -2 $DIR_RAW/${ID_FATHER}_2.fastq.gz \
		 -E $EMAIL \
		 -G GSC -s $ID_FATHER --readlength $READ_LENGTH -S PBS \
		 --mtoolbox $MTOOLBOX_CONFIG_FILE --metrics-exome


## STEP 2) processes and annotate vcfs and convert to GEMINI database, creates "Bam2Gemini" script
python $SCRIPT_BAM2GEMINI -G GSC -T Exome \
        -v GVCF -V ${ID_PROBAND}_dupremoved_realigned_HaplotypeCaller.g.vcf,${ID_MOTHER}_dupremoved_realigned_HaplotypeCaller.g.vcf,${ID_FATHER}_dupremoved_realigned_HaplotypeCaller.g.vcf \
        -P $ID_FAMILY.ped -p 8 -m 30G -F $ID_FAMILY \
 		-E $EMAIL \
	-D $DB_DIR \
	-A $BASE_DIR \
        -d $DIR_WORKING


## STEP 3) make a PED file
python $PEDMAKER \
        --proband ${ID_PROBAND},${PROBAND_SEX},affected \
        --mother ${ID_MOTHER},female,unaffected \
        --father ${ID_FATHER},male,unaffected \
        --FamilyID $ID_FAMILY \
        -O $DIR_WORKING/$ID_FAMILY.ped

