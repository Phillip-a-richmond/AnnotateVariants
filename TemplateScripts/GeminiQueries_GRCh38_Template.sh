#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

## CPU Usage
## 16 Gb of RAM for the whole job
#SBATCH --mem=16G

## Using 2 CPUs
#SBATCH --cpus-per-task=2

## Running for a max time of 8 hours
#SBATCH --time=8:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


# Gemini Query script
# Phillip Richmond
# ChangeLog:
# July 10th 2018, script created with addition of gnomad exome and revised filters
# August 30th 2018, added PrimateAI to query, fixed ClinVar to be split to pathogenic/likely pathogenic, fixed General Damaging to be split with Het/Homo
# September 10th, 2018 - Fixed issue with general damaging which impacted reporting of variants
# NOTES:
# Updates to this script should accomodate differences in how many individuals are present
# March 21st, 2019 - Major update:
	# Removed LOOSE DP/GQ thresholds
	# Added SPLICEAI
	# Redid columns
	# Added FATHMM-XF-NONCODING
	# Added NONCODING - not coding, CADD >= 20 or FATHMM-XF >= 0.9
	# Added pop threshold on Clinvar to get rid of actionable polymorphisms
	# Added FATHMM-XF-NONCODING
	# Added De Novo LOW (de novo variants without predicted HIGH/MED impact, useful for WGS
	# Added new GEMINI location explicitly
	# Added mendel_errors
	# Changed strict min DP to 15
# May 23, 2019:
	# Set maximum number of homozygotes for cph query to 15
# Overhaul for GPCC cluster, GRCh38 version
# September 11, 2020 (never forget)
	# GRCh38: CADDv1.4, CADDv1.6, gnomADv3, clinvar, 
# March 9th, 2021: Update for GPCC
	# functional update for GPCC, including pairing with an update to the tableannotator 

# Load the environment
ANNOTATE_VARIANTS_DIR=annotate_variants_dir
source $ANNOTATE_VARIANTS_DIR/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATE_VARIANTS_DIR/opt/AnnotateVariantsEnvironment/

## Change to working directory
FAMILY_ID=family_id
GEMINIDB=${FAMILY_ID}_Gemini.db
WORKING_DIR=working_dir
DATA_DIR=$WORKING_DIR
cd $WORKING_DIR

# Here I'm pulling with singularity, and recoding the variable GEMINI to equal the call to Gemini in the singularity container
# This assumes a data directory has been specified, as seen above
module load singularity
export SINGULARITY_CACHEDIR=$PWD
singularity pull docker://quay.io/biocontainers/gemini:0.20.0--py27_0
GEMINI="singularity run -B /usr/lib/locale/:/usr/lib/locale/ -B $DATA_DIR:/data/ gemini_0.20.0--py27_0.sif  gemini" 
MTOOLBOX_RSCRIPT=$ANNOTATE_VARIANTS_DIR/MToolBox_config_files/Mtoolbox.R
TableAnnotator=$ANNOTATE_VARIANTS_DIR/TableAnnotators/GeminiTable2CVL.py

# Define Inheritance Model Variant Files
# Phase 1
AUTODOM_OUT=$WORKING_DIR${FAMILY_ID}_autoDom
DENOVO_OUT=$WORKING_DIR${FAMILY_ID}_deNovo
DENOVO_LOW_OUT=$WORKING_DIR${FAMILY_ID}_deNovoLow
RECESSIVE_OUT=$WORKING_DIR${FAMILY_ID}_recessive
COMPOUND_HET_OUT=$WORKING_DIR${FAMILY_ID}_compoundHet
X_RECESSIVE_OUT=$WORKING_DIR${FAMILY_ID}_Xrecessive
X_DOMINANT_OUT=$WORKING_DIR${FAMILY_ID}_Xdominant
X_DENOVO_OUT=$WORKING_DIR${FAMILY_ID}_Xdenovo
MENDEL_ERRORS_OUT=$WORKING_DIR${FAMILY_ID}_mendelErrors

# Phase 2
DBSTATFILE=${FAMILY_ID}_GeminiDB_Stats.txt

# Phase 3
GENERAL_DAMAGING_HET=${FAMILY_ID}_GeneralDamaging_Het
GENERAL_DAMAGING_HOMO=${FAMILY_ID}_GeneralDamaging_Homo
CLINVAR_HITS=${FAMILY_ID}_Clinvar_Hits
SPLICING_HITS=${FAMILY_ID}_SpliceCandidates
NONCODING_HITS=${FAMILY_ID}_NoncodingCandidates

# Columns to present within CVL
COLUMNS="chrom, start, end, ref, alt, gene, clinvar_variation_id, ensembl_gene_id, exon, aa_change, codon_change, vepprotein_position, transcript, biotype, hgvsc, hgvsp, impact, impact_severity, rs_ids, filter, gts, gt_ref_depths, gt_alt_depths, gt_alt_freqs, gt_quals, gt_types, gt_phases, gt_depths, gt_quals, gt_alt_freqs, gnomad_genome_ac_global, gnomad_genome_an_global, gnomad_genome_af_global, gnomad_genome_hom_global, gnomad_popmax_af, gnomad_nhomalt, cadd_v1_4, cadd_v1_4_indel, cadd_v1_6, cadd_v1_6_indel, polyphen_pred, polyphen_score, sift_pred, sift_score, vepspliceai_pred_dp_ag, vepspliceai_pred_dp_al, vepspliceai_pred_dp_dg, vepspliceai_pred_dp_dl, vepspliceai_pred_ds_ag, vepspliceai_pred_ds_al, vepspliceai_pred_ds_dg, vepspliceai_pred_ds_dl, vepspliceai_pred_symbol, vepmaxentscan_alt, vepmaxentscan_diff, vepmaxentscan_ref, veppubmed, vepexisting_variation, vepclin_sig, clinvar_dbinfo, clinvar_disease_name, clinvar_pathogenic,cosmic_coding_ids, cosmic_count_observed, somatic, clin_sig, pheno,hgvs_offset,impact_so"

# Some variables for filters
CODING='is_coding=1'
EXONIC='is_exonic=1'
SPLICING='is_splicing=1'

LOF='is_lof=1'

IMPACT_HIGH="impact_severity=='HIGH'"
IMPACT_MED="impact_severity=='MED'"
IMPACT_LOW="impact_severity=='LOW'"

FILTER='(filter is NULL)'

CADD='((cadd_v1_4 >= 20) OR (cadd_v1_4_indel >= 20) OR (cadd_v1_6 >= 20) OR (cadd_v1_6_indel >= 20)) '
NONCODING="($CADD) AND (NOT $CODING)"
SPLICEAI='( (vepspliceai_pred_ds_ag > 0.2) OR ( vepspliceai_pred_ds_al > 0.2 ) OR ( vepspliceai_pred_ds_dl > 0.2 ) OR ( vepspliceai_pred_ds_dg > 0.2 ) )'
STRICT_MIN_DP=15
STRICT_MIN_GQ=30


#########################################################################################

#####################
#       Phase 1     #
#####################
echo "Starting Phase 1"
# Inheritance models
## Recessive
$GEMINI autosomal_recessive \
	--columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)" \
	$GEMINIDB > $RECESSIVE_OUT
	
python $TableAnnotator -i $RECESSIVE_OUT -o ${RECESSIVE_OUT}_annotated.txt


# Compound Het Variants, number of homozygotes == 15
$GEMINI comp_hets \
	--columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)" \
	$GEMINIDB > $COMPOUND_HET_OUT
python $TableAnnotator -i $COMPOUND_HET_OUT -o ${COMPOUND_HET_OUT}_annotated.txt

# De novo
$GEMINI de_novo \
	--columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)" \
	$GEMINIDB > $DENOVO_OUT
python $TableAnnotator -i $DENOVO_OUT -o ${DENOVO_OUT}_annotated.txt

# De novo Not HIGH/MED
$GEMINI de_novo \
	--columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_LOW)" \
	$GEMINIDB > $DENOVO_LOW_OUT
python $TableAnnotator -i $DENOVO_LOW_OUT -o ${DENOVO_LOW_OUT}_annotated.txt

# X Dominant
$GEMINI x_linked_dominant  \
	--columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	$GEMINIDB > $X_DOMINANT_OUT
python $TableAnnotator -i $X_DOMINANT_OUT -o ${X_DOMINANT_OUT}_annotated.txt

# X De Novo
$GEMINI x_linked_de_novo \
	--columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	$GEMINIDB > $X_DENOVO_OUT
python $TableAnnotator -i $X_DENOVO_OUT -o ${X_DENOVO_OUT}_annotated.txt

# X Recessive
$GEMINI x_linked_recessive \
        --columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)"\
        $GEMINIDB > $X_RECESSIVE_OUT
python $TableAnnotator -i $X_RECESSIVE_OUT -o ${X_RECESSIVE_OUT}_annotated.txt

# Autosomal Dominant
$GEMINI autosomal_dominant \
	--columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	$GEMINIDB > $AUTODOM_OUT
python $TableAnnotator -i $AUTODOM_OUT -o ${AUTODOM_OUT}_annotated.txt

# Mendel Errors
$GEMINI mendel_errors \
	--columns "$COLUMNS" \
	--filter "$FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	$GEMINIDB > $MENDEL_ERRORS_OUT
python $TableAnnotator -i $MENDEL_ERRORS_OUT -o ${MENDEL_ERRORS_OUT}_annotated.txt

##########################################################################################


#####################
#       Phase 2     #
#####################

# General Queries to the database:
GNOMAD_GENOME_COMMON='(gnomad_genome_af_global > 0.01)'
GNOMAD_GENOME_RARE='((gnomad_genome_af_global <= 0.01) or (gnomad_genome_af_global is NULL)) '

DBSNP='rs_ids is not NULL'

UTR='( (is_coding = 0) and (is_exonic = 1))'
# Check your sample info
$GEMINI query -q "SELECT * FROM samples" $GEMINIDB > $DBSTATFILE


#I will print out a table
echo "Annotation	Total	Common	Rare" >> $DBSTATFILE

# Total number of variants in the file
TOTAL_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants WHERE $FILTER" $GEMINIDB`
rareTOTAL_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants WHERE $FILTER AND $GNOMAD_GENOME_RARE" $GEMINIDB`
commonTOTAL_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants WHERE $FILTER AND ($GNOMAD_GENOME_COMMON)" $GEMINIDB`


echo "TOTAL	$TOTAL_VARS	$commonTOTAL_VARS	$rareTOTAL_VARS" >> $DBSTATFILE

# UTR vars
UTR_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $UTR AND $FILTER" $GEMINIDB`
rareUTR_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $UTR AND $FILTER AND $GNOMAD_GENOME_RARE "  $GEMINIDB`
commonUTR_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $UTR AND $FILTER AND ($GNOMAD_GENOME_COMMON)"  $GEMINIDB`


echo "UTR	$UTR_VARS	$commonUTR_VARS	$rareUTR_VARS" >> $DBSTATFILE

# Exonic variants
EXONIC_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $FILTER" $GEMINIDB`
rareEXONIC_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $FILTER AND $GNOMAD_GENOME_RARE"  $GEMINIDB`
commonEXONIC_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $FILTER AND $GNOMAD_GENOME_COMMON"  $GEMINIDB`

echo "EXONIC	$EXONIC_VARS	$commonEXONIC_VARS	$rareEXONIC_VARS" >> $DBSTATFILE

# Coding Variants
CODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $FILTER" $GEMINIDB`
rareCODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $FILTER AND $GNOMAD_GENOME_RARE"  $GEMINIDB`
commonCODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $FILTER AND ($GNOMAD_GENOME_COMMON)"  $GEMINIDB`

echo "CODING	$CODING_VARS	$commonCODING_VARS	$rareCODING_VARS" >> $DBSTATFILE

#Splicing variants
SPLICING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $FILTER" $GEMINIDB`
rareSPLICING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $FILTER AND $GNOMAD_GENOME_RARE"  $GEMINIDB`
commonSPLICING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $FILTER AND ($GNOMAD_GENOME_COMMON)"  $GEMINIDB`

echo "SPLICING	$SPLICING_VARS	$commonSPLICING_VARS	$rareSPLICING_VARS" >> $DBSTATFILE

# LOF variants
LOF_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $FILTER" $GEMINIDB`
rareLOF_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $FILTER AND $GNOMAD_GENOME_RARE"  $GEMINIDB`
commonLOF_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $FILTER AND ($GNOMAD_GENOME_COMMON)"  $GEMINIDB`

echo "LOF	$LOF_VARS	$commonLOF_VARS	$rareLOF_VARS" >> $DBSTATFILE

# CADD >= 20
CADD_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $CADD AND $FILTER" $GEMINIDB`
rareCADD_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $CADD AND $FILTER AND $GNOMAD_GENOME_RARE"  $GEMINIDB`
commonCADD_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $CADD AND $FILTER AND ($GNOMAD_GENOME_COMMON)"  $GEMINIDB`

echo "CADD	$CADD_VARS	$commonCADD_VARS	$rareCADD_VARS" >> $DBSTATFILE



#####################
#       Phase 3     #
#####################
#GNOMAD_GENOME_RARE="( (gnomad_genome_af_global <= 0.01 or gnomad_genome_af_global is NULL) AND (gnomad_genome_hom_global <= 10 or gnomad_genome_hom_global is NULL) )"
#GNOMAD_EXOME_RARE="( (gnomad_exome_af_global <= 0.01 or gnomad_exome_af_global is NULL) AND (gnomad_exome_hom_global <= 10 or gnomad_exome_hom_global is NULL) )"
# NOTE: Changing this allows homozygotes to escape into the general damaging. we are fine with that.
GNOMAD_GENOME_RARE='((gnomad_genome_af_global <= 0.01) or (gnomad_genome_af_global is NULL)) '

# General Damaging Variant List for the proband, which will include clinvar pathogenic hits no matter what
# note: Added rarity to clinvar here
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE ( (clinvar_pathogenic == 'Pathogenic') OR (clinvar_pathogenic == 'Likely Pathogenic') ) AND (clinvar_disease_name is not NULL) AND $FILTER AND $GNOMAD_GENOME_RARE" \
	--gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(any)" \
	--header \
	$GEMINIDB > $CLINVAR_HITS

python  $TableAnnotator -i $CLINVAR_HITS -o ${CLINVAR_HITS}_annotated.txt

# General damaging het (HIGH/MED impact)
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	--gt-filter "(gt_types).(phenotype==2).(==HET).(any)" \
	--header \
	$GEMINIDB > $GENERAL_DAMAGING_HET

python $TableAnnotator -i $GENERAL_DAMAGING_HET -o ${GENERAL_DAMAGING_HET}_annotated.txt

# General damaging homo (HIGH/MED impact)
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER" \
	--gt-filter "(gt_types).(phenotype==2).(==HOM_ALT).(any)" \
	--header \
	$GEMINIDB > $GENERAL_DAMAGING_HOMO

python $TableAnnotator -i $GENERAL_DAMAGING_HOMO -o ${GENERAL_DAMAGING_HOMO}_annotated.txt

# Possible Intronic
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE $SPLICEAI AND $FILTER" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(any)" \
        --header \
        $GEMINIDB > $SPLICING_HITS
python $TableAnnotator -i $SPLICING_HITS -o ${SPLICING_HITS}_annotated.txt

# Noncoding high impact
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE $NONCODING AND $FILTER" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(any)" \
        --header \
        $GEMINIDB > $NONCODING_HITS
python $TableAnnotator -i $NONCODING_HITS -o ${NONCODING_HITS}_annotated.txt

#####################
#       Phase 4     #
#####################


# Integrate into single file
# Header
CVL=$WORKING_DIR${FAMILY_ID}_CVL.tsv
DATE=`date`
echo "FAMILY:	$FAMILY_ID" > $CVL
echo "DATE CVL PRODUCED:	$DATE" >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# General Stats
cat $DBSTATFILE >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# De novo
echo "De Novo (HIGH/MED)" >> $CVL
cat ${DENOVO_OUT}_annotated.txt >> $CVL

# De novo LOW
echo "De Novo (LOW)" >> $CVL
cat ${DENOVO_OUT_LOW}_annotated.txt >> $CVL

echo "---------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Recessive
echo "Recessive" >> $CVL
cat ${RECESSIVE_OUT}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Compound Het
echo "Compound Het" >> $CVL
cat ${COMPOUND_HET_OUT}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Dominant
echo "Autosomal Dominant" >> $CVL
cat ${AUTODOM_OUT}_annotated.txt >> $CVL


echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# X de novo
echo "X De Novo" >> $CVL
cat ${X_DENOVO_OUT}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL


#X Recessive
echo "X Recessive" >> $CVL
cat ${X_RECESSIVE_OUT}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL


#X Autosomal dominant
echo "X Autosomal Dominant" >> $CVL
cat ${X_DOMINANT_OUT}_annotated.txt >> $CVL


echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Mendelian Errors
echo "Mendelian Errors" >> $CVL
cat ${MENDEL_ERRORS_OUT}_annotated.txt >> $CVL
echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

############
# MToolBox #
############

cp $MTOOLBOX_RSCRIPT ./
Rscript ./Mtoolbox.R
echo "MToolbox Output" >> $CVL
cat MToolBox_annotated.txt >> $CVL
echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL

# General Damaging
echo "General Damaging Homozygous" >> $CVL
cat ${GENERAL_DAMAGING_HOMO}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

echo "General Damaging Heterozygous" >> $CVL
cat ${GENERAL_DAMAGING_HET}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Clinvar Hits
echo "Clinvar Hits" >> $CVL
cat ${CLINVAR_HITS}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

#######################
# Other, deeper dives #
#######################

# Splicing Candidates

echo "Splicing Candidates">> $CVL
cat ${SPLICING_HITS}_annotated.txt >> $CVL
echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

## Noncoding Others with CADD / FATHMM-XF
#echo "Noncoding Candidates, cis-regulatory" >> $CVL
#cat ${NONCODING_HITS}_annotated.txt >> $CVL
#echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
############

#########################
# Phase 5 bis for macro #
#########################


# Integrate into single file that can be run in excel macro to create tabs
# Header
CVL_macro=$WORKING_DIR${FAMILY_ID}_CVL_macro.tsv
DATE=`date`
echo "Nothing	FAMILY:	$FAMILY_ID" > $CVL_macro
echo "Nothing	DATE CVL PRODUCED:	$DATE" >> $CVL_macro
echo "Nothing	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# General Stats
#cat $DBSTATFILE >> $CVL_macro
perl -ne 'print "Nothing\t$_" unless /^gene\s+Gene_Name\s+HPO/' $DBSTATFILE >> $CVL_macro

echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# De novo
echo "De novo	De novo" >> $CVL_macro
#cat ${DENOVO_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "De novo\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${DENOVO_OUT}_annotated.txt >> $CVL_macro

echo "De novo" >> $CVL_macro
echo "De novo	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL_macro
echo "De novo" >> $CVL_macro

# De novo Low
echo "De novo	De novo low" >> $CVL_macro
#cat ${DENOVO_LOW_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "De novo\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${DENOVO_LOW_OUT}_annotated.txt >> $CVL_macro


echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# Recessive
echo "Autosomal Rec (Homozygous)	Autosomal Rec (Homozygous)" >> $CVL_macro
#cat ${RECESSIVE_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "Autosomal Rec (Homozygous)\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${RECESSIVE_OUT}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# Compound Het
echo "Recessive (Compound Het)	Compound Het" >> $CVL_macro
#cat ${COMPOUND_HET_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "Recessive (Compound Het)\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${COMPOUND_HET_OUT}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# Dominant
echo "Autosomal Dom (Heterozygous)	Autosomal Dominant" >> $CVL_macro
#cat ${AUTODOM_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "Autosomal Dom (Heterozygous)\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${AUTODOM_OUT}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# X de novo
echo "X De Novo	X De Novo" >> $CVL_macro
#cat ${X_DENOVO_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "X De Novo\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${X_DENOVO_OUT}_annotated.txt >> $CVL_macro

echo "X De Novo" >> $CVL_macro
echo "X De Novo	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL_macro
echo "X De Novo" >> $CVL_macro

echo "X De Novo	X De Novo loose" >> $CVL_macro
#cat ${X_DENOVO_OUT_LOOSE}_annotated.txt >> $CVL_macro
perl -ne 'print "X De Novo\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${X_DENOVO_OUT_LOOSE}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

#X Recessive
echo "X Recessive	X Recessive" >> $CVL_macro
#cat ${X_RECESSIVE_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "X Recessive\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${X_RECESSIVE_OUT}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

#X Autosomal dominant
echo "X Dominant	X Autosomal Dominant" >> $CVL_macro
#cat ${X_DOMINANT_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "X Dominant\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${X_DOMINANT_OUT}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# Mendelian Errors
echo "Mendelian errors	Mendelian Errors" >> $CVL_macro
#cat ${MENDEL_ERRORS_OUT}_annotated.txt >> $CVL_macro
perl -ne 'print "Mendelian errors\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${MENDEL_ERRORS_OUT}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

############
# MToolBox #
############

cp $MTOOLBOX_RSCRIPT ./
Rscript ./Mtoolbox.R
echo "Mitochondrial	MToolbox Output" >> $CVL_macro
#cat MToolBox_annotated.txt >> $CVL_macro
perl -ne 'print "Mitochondrial\t$_" unless /^.Variant.Allele\s+Samples\s+HF/' MToolBox_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# General Damaging
echo "Damaging (No Inheritance)	General Damaging Homozygous" >> $CVL_macro
#cat ${GENERAL_DAMAGING_HOMO}_annotated.txt >> $CVL_macro
perl -ne 'print "Damaging (No Inheritance)\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${GENERAL_DAMAGING_HOMO}_annotated.txt >> $CVL_macro

echo "Damaging (No Inheritance)" >> $CVL_macro
echo "Damaging (No Inheritance)	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL_macro
echo "Damaging (No Inheritance)" >> $CVL_macro

echo "Damaging (No Inheritance)	General Damaging Heterozygous" >> $CVL_macro
#cat ${GENERAL_DAMAGING_HET}_annotated.txt >> $CVL_macro
perl -ne 'print "Damaging (No Inheritance)\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${GENERAL_DAMAGING_HET}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# Clinvar Hits
echo "ClinVar Hits (No inheritance)	Clinvar Hits" >> $CVL_macro
#cat ${CLINVAR_HITS}_annotated.txt >> $CVL_macro
perl -ne 'print "ClinVar Hits (No inheritance)\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${CLINVAR_HITS}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

#######################
# Other, deeper dives #
#######################

# Splicing Candidates

echo "Splicing Candidates	Splicing Candidates">> $CVL_macro
#cat ${SPLICING_HITS}_annotated.txt >> $CVL_macro
perl -ne 'print "Splicing Candidates\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${SPLICING_HITS}_annotated.txt >> $CVL_macro

echo "Nothing" >> $CVL_macro
echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
echo "Nothing" >> $CVL_macro

# Noncoding Others with CADD / FATHMM-XF
#echo "Non Coding - Cis Regulatory	Noncoding Candidates, cis-regulatory" >> $CVL_macro
#cat ${NONCODING_HITS}_annotated.txt >> $CVL_macro
#perl -ne 'print "Non Coding - Cis Regulatory\t$_" unless /^gene\s+Gene_Name\s+HPO/' ${NONCODING_HITS}_annotated.txt >> $CVL_macro

#echo "Nothing" >> $CVL_macro
#echo "Nothing	==========================================================================================================================================================================================================================================" >> $CVL_macro
#echo "Nothing" >> $CVL_macro

############
# Clean Up
rm *_annotated.txt
rm $AUTODOM_OUT
rm $DENOVO_OUT
rm $DENOVO_LOW_OUT
rm $RECESSIVE_OUT
rm $COMPOUND_HET_OUT
rm $X_RECESSIVE_OUT
rm $X_DOMINANT_OUT
rm $X_DENOVO_OUT
rm $AUTODOM_OUT_LOOSE
rm $DENOVO_OUT_LOOSE
rm $RECESSIVE_OUT_LOOSE
rm $COMPOUND_HET_OUT_LOOSE
rm $X_RECESSIVE_OUT_LOOSE
rm $X_DOMINANT_OUT_LOOSE
rm $X_DENOVO_OUT_LOOSE
rm $GENERAL_DAMAGING_HET
rm $GENERAL_DAMAGING_HOMO
rm $CLINVAR_HITS
rm $SPLICING_HITS
rm $NONCODING_HITS
rm $MENDEL_ERRORS_OUT
