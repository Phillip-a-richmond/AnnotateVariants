#!/bin/bash

#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########

# Load environment (need this for outlier with ehdn)
MINICONDA_DIR=miniconda_dir
source $MINICONDA_DIR/opt/miniconda2/etc/profile.d/conda.sh
conda activate EHdn


# Number of threads
NSLOTS=$SLURM_CPUS_PER_TASK

# Go to the submission directory (where the sbatch was entered)
cd $SLURM_SUBMIT_DIR
WORKING_DIR=working_dir

## Set working space
mkdir -p $WORKING_DIR
cd $WORKING_DIR

#### GRCh38 #### 
echo "GRCh38 genome"
GENOME=genome_build
FASTA_DIR=fasta_dir
FASTA_FILE=fasta_file

BAM_DIR=$WORKING_DIR
FAMILY_ID=family_id
PROBAND_SAMPLEID=proband_id
MOTHER_SAMPLEID=mother_id
FATHER_SAMPLEID=father_id

PROBAND_BAM=${PROBAND_SAMPLEID}_$GENOME.dupremoved.sorted.bam
FATHER_BAM=${FATHER_SAMPLEID}_$GENOME.dupremoved.sorted.bam
MOTHER_BAM=${MOTHER_SAMPLEID}_$GENOME.dupremoved.sorted.bam

PROBAND_BASE=${PROBAND_SAMPLEID}_EHdn
FATHER_BASE=${FATHER_SAMPLEID}_EHdn
MOTHER_BASE=${MOTHER_SAMPLEID}_EHdn


# Run ExpansionHunter DeNovo 
EHDN_DIR=ehdn_dir_var
EHDN=${EHDN_DIR}/bin/ExpansionHunterDenovo
ANNOTATE_EHDN=${EHDN_DIR}/scripts/annotate_ehdn.sh
ANNOVAR_DIR=annovar_dir
ANNOVAR_ANNOTATE_VARIATION=$ANNOVAR_DIR/annotate_variation.pl
ANNOVAR_HUMANDB=$ANNOVAR_DIR/humandb/
ANNOVAR_GENOMEVERSION=hg38
BACKGROUND_DIR=ehdn_background_dir


# Step 1 - Run EHDN profile
# For now, we will ignore the parental calls
$EHDN profile \
	--reads $PROBAND_BAM \
	--reference $FASTA_DIR/$FASTA_FILE \
	--output-prefix $PROBAND_BASE

#$EHDN profile \
#	--reads $FATHER_BAM \
#	--reference $FASTA_DIR/$FASTA_FILE \
#	--output-prefix $FATHER_BASE
#
#
#$EHDN profile \
#	--reads $MOTHER_BAM \
#	--reference $FASTA_DIR/$FASTA_FILE \
#	--output-prefix $MOTHER_BASE

# Step 2: Combine counts, perform outlier tests, annotate
JSON_FILES=( $WORKING_DIR/*.str_profile.json ) 
for SAMPLE_JSON in ${JSON_FILES[@]}
	do
	# A) Make Manifest
	IFS='/' read -a array <<< $SAMPLE_JSON
	SJSON=${array[-1]}
	IFS='.' read -a array2 <<< $SJSON
	SAMPLE_ID=${array2[0]}
	IFS='_' read -a array2 <<< $SAMPLE_ID
	MANIFEST=${WORKING_DIR}${SAMPLE_ID}_manifest.tsv
	rm $MANIFEST
	
	# Used for re-running if failed
#	if [ -f ${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Locus.txt ]; then
#	        exit
#	fi
	
	#echo $SJSON
#	echo $SAMPLE_ID
	echo "$SAMPLE_ID	case	${SAMPLE_JSON}" >> $MANIFEST
	for BACKGROUND_JSON in $(ls $BACKGROUND_DIR/*json)
	do
		#echo $BACKGROUND_JSON
		# Split the full filepath, only take last file (ignoring directory)
		IFS='/' read -a array1 <<< $BACKGROUND_JSON
		BJSON=${array1[-1]}
		# Split the filename, removing JSON, could also substring but whatever
		IFS='.' read -a array2 <<< $BJSON
		BJSON_SAMPLE_ID=${array2[0]}
		if ! [  $SAMPLE_ID == $BJSON_SAMPLE_ID ];
		then
#			echo $BJSON_SAMPLE_ID
#			echo $SAMPLE_ID
			echo "$BJSON_SAMPLE_ID	control	$BACKGROUND_JSON" >> $MANIFEST
		fi
	done
	
	
	# B) Combine counts
	$EHDN merge --output-prefix $SAMPLE_ID \
	--reference $FASTA_DIR/$FASTA_FILE \
	--manifest $MANIFEST 	
	# For reference, this creates a file with combined counts called: ${SAMPLE_ID}.multisample_profile.json


	# Step 3A: Run Outlier Analysis at the locus level
	OUTPUT_RESULTS_LOCUS=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Locus.txt
	$PYTHON $EHDN_DIR/scripts/outlier.py locus  \
	--manifest $MANIFEST \
	--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
	--output $OUTPUT_RESULTS_LOCUS
	
	# And Sort the results
	(head -n1 $OUTPUT_RESULTS_LOCUS && tail -n +2 $OUTPUT_RESULTS_LOCUS | sort -k5Vr ) > ${OUTPUT_RESULTS_LOCUS}.sorted.tsv
	mv ${OUTPUT_RESULTS_LOCUS}.sorted.tsv $OUTPUT_RESULTS_LOCUS
	
	# Step 3B: Run Outlier Analysis at the motif level
	OUTPUT_RESULTS_MOTIF=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Motif.txt
	$PYTHON $EHDN_DIR/scripts/outlier.py motif  \
	--manifest $MANIFEST \
	--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
	--output $OUTPUT_RESULTS_MOTIF
	
	# And Sort the results
	(head -n1 $OUTPUT_RESULTS_MOTIF && tail -n +2 $OUTPUT_RESULTS_MOTIF | sort -k2Vr ) > ${OUTPUT_RESULTS_MOTIF}.sorted.tsv
	mv ${OUTPUT_RESULTS_MOTIF}.sorted.tsv $OUTPUT_RESULTS_MOTIF
	
	
	# Step 4, Annotate with ANNOVAR
	
	sh $ANNOTATE_EHDN \
		--ehdn-results $OUTPUT_RESULTS_LOCUS \
		--ehdn-annotated-results ${OUTPUT_RESULTS_LOCUS::-4}_Annotated.tsv \
		--annovar-annotate-variation $ANNOVAR_ANNOTATE_VARIATION \
		--annovar-humandb $ANNOVAR_HUMANDB \
		--annovar-buildver $ANNOVAR_GENOMEVERSION 
done	

