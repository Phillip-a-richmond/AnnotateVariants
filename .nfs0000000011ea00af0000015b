# Gemini Query script
# Phillip Richmond
# ChangeLog:
# July 10th 2018, script created with addition of gnomad exome and revised filters

# NOTES: 
# Updates to this script should accomodate differences in how many individuals are present

## Define variables
WORKING_DIR='/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/T242/'
GEMINIDB='T242.db'
FAMILY_ID='T242'
TableAnnotator=/mnt/causes-vnx1/Pipelines/AnnotateVariants/TableAnnotators/GeminiTable2CVL.py

source /opt/tools/hpcenv.sh



## Change to working directory
cd $WORKING_DIR

# Define Inheritance Model Variant Files
# Phase 1
AUTODOM_OUT=$WORKING_DIR${FAMILY_ID}_autoDom
DENOVO_OUT=$WORKING_DIR${FAMILY_ID}_deNovo
RECESSIVE_OUT=$WORKING_DIR${FAMILY_ID}_recessive
COMPOUND_HET_OUT=$WORKING_DIR${FAMILY_ID}_compoundHet
X_RECESSIVE_OUT=$WORKING_DIR${FAMILY_ID}_Xrecessive
X_DOMINANT_OUT=$WORKING_DIR${FAMILY_ID}_Xdominant
X_DENOVO_OUT=$WORKING_DIR${FAMILY_ID}_Xdenovo

# Phase 2
AUTODOM_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_autoDom_loose
DENOVO_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_deNovo_loose
RECESSIVE_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_recessive_loose
COMPOUND_HET_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_compoundHet_loose
X_RECESSIVE_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_Xrecessive_loose
X_DOMINANT_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_Xdominant_loose
X_DENOVO_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_Xdenovo_loose

# Phase 3
DBSTATFILE=${FAMILY_ID}_GeminiDB_Stats.txt

# Phase 4
GENERAL_DAMAGING=${FAMILY_ID}_GeneralDamaging
GENERAL_DAMAGING_LOOSE=${FAMILY_ID}_GeneralDamaging_loose
CLINVAR_HITS=${FAMILY_ID}_Clinvar_Hits

# Columns to present within CVL
COLUMNS="chrom, start, end, ref, alt, gene, exon, aa_change, impact, impact_severity, rs_ids, filter, gts, gt_ref_depths, gt_alt_depths, gt_quals, in_segdup, confidentregion, inhousedb_ac, aaf_gnomad_exome_all, gnomad_exome_hom_all, gnomad_exome_af_afr, gnomad_exome_hom_afr, gnomad_exome_af_amr, gnomad_exome_hom_amr, gnomad_exome_af_asj, gnomad_exome_hom_asj, gnomad_exome_af_eas, gnomad_exome_hom_eas, gnomad_exome_af_fin, gnomad_exome_hom_fin, gnomad_exome_af_nfe, gnomad_exome_hom_nfe, gnomad_exome_af_oth, gnomad_exome_hom_oth, aaf_gnomad_genome_all, gnomad_genome_hom_all, gnomad_genome_af_afr, gnomad_genome_hom_afr, gnomad_genome_af_amr, gnomad_genome_hom_amr, gnomad_genome_af_asj, gnomad_genome_hom_asj, gnomad_genome_af_eas, gnomad_genome_hom_eas, gnomad_genome_af_fin, gnomad_genome_hom_fin, gnomad_genome_af_nfe, gnomad_genome_hom_nfe, gnomad_genome_af_oth, gnomad_genome_hom_oth, cadd, pp2hdiv, pp2hvar, clinvar_disease_name, clinvar_pathogenic, clinvar_dbinfo"

# Define Variant Annotation Cutoffs
GNOMAD_GENOME_RARE='( (aaf_gnomad_genome_all <= 0.01 or aaf_gnomad_genome_all is NULL) AND (gnomad_genome_hom_all <= 10 or gnomad_genome_hom_all is NULL) )'
GNOMAD_EXOME_RARE='( (aaf_gnomad_exome_all <= 0.01 or aaf_gnomad_exome_all is NULL) AND (gnomad_exome_hom_all <= 10 or gnomad_exome_hom_all is NULL) )'
INHOUSE_RARE='(inhousedb_ac <= 3 or inhousedb_ac is NULL)'

CODING='is_coding=1'
EXONIC='is_exonic=1'
SPLICING='is_splicing=1'

LOF='is_lof=1'

IMPACT_HIGH="impact_severity=='HIGH'"
IMPACT_MED="impact_severity=='MED'"
IMPACT_LOW="impact_severity=='LOW'"

FILTER='filter is NULL'
CONFIDENTREGION='confidentregion = 1'
SEGDUP='in_segdup=0'

STRICT_MIN_DP=20
STRICT_MIN_GQ=30

LOOSE_MIN_DP=15
LOOSE_MIN_GQ=20

# STRICT FILTER
# --filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\

# LOOSE FILTER
# --filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"

##########################################################################################

#####################
#       Phase 1     #
#####################
# Inheritance models
## Recessive
gemini autosomal_recessive \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $RECESSIVE_OUT
python $TableAnnotator -i $RECESSIVE_OUT -o ${RECESSIVE_OUT}_annotated.txt

# Compound Het Variants
gemini comp_hets \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $COMPOUND_HET_OUT
python $TableAnnotator -i $COMPOUND_HET_OUT -o ${COMPOUND_HET_OUT}_annotated.txt

# De novo
gemini de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $DENOVO_OUT
python $TableAnnotator -i $DENOVO_OUT -o ${DENOVO_OUT}_annotated.txt


# X Dominant
gemini x_linked_dominant  \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $X_DOMINANT_OUT
python $TableAnnotator -i $X_DOMINANT_OUT -o ${X_DOMINANT_OUT}_annotated.txt

# X De Novo
gemini x_linked_de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $X_DENOVO_OUT
python $TableAnnotator -i $X_DENOVO_OUT -o ${X_DENOVO_OUT}_annotated.txt

# X Recessive
gemini x_linked_recessive \
        --columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
        -d $STRICT_MIN_DP \
        --min-gq $STRICT_MIN_GQ \
        $GEMINIDB > $X_RECESSIVE_OUT
python $TableAnnotator -i $X_RECESSIVE_OUT -o ${X_RECESSIVE_OUT}_annotated.txt

# Autosomal Dominant
gemini autosomal_dominant \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $LOOSE_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $AUTODOM_OUT
python $TableAnnotator -i $AUTODOM_OUT -o ${AUTODOM_OUT}_annotated.txt


##########################################################################################

#####################
#       Phase 2     #
#####################
# Repeat but with Looser filters, specifically on genotype quality, depth, and whether or not it's in a segdup

# Inheritance models
## Recessive
gemini autosomal_recessive \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	-d $LOOSE_MIN_DP \
	--min-gq $LOOSE_MIN_GQ \
	$GEMINIDB > $RECESSIVE_OUT_LOOSE
python $TableAnnotator -i $RECESSIVE_OUT_LOOSE -o ${RECESSIVE_OUT_LOOSE}_annotated.txt


# Compound Het Variants
gemini comp_hets \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	-d $LOOSE_MIN_DP \
	--min-gq $LOOSE_MIN_GQ \
	$GEMINIDB > $COMPOUND_HET_OUT_LOOSE
python $TableAnnotator -i $COMPOUND_HET_OUT_LOOSE -o ${COMPOUND_HET_OUT_LOOSE}_annotated.txt

# De novo
gemini de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	-d $LOOSE_MIN_DP \
	--min-gq $LOOSE_MIN_GQ \
	$GEMINIDB > $DENOVO_OUT_LOOSE
python $TableAnnotator -i $DENOVO_OUT_LOOSE -o ${DENOVO_OUT_LOOSE}_annotated.txt

# X Dominant
gemini x_linked_dominant  \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	-d $LOOSE_MIN_DP \
	--min-gq $LOOSE_MIN_GQ \
	$GEMINIDB > $X_DOMINANT_OUT_LOOSE
python $TableAnnotator -i $X_DOMINANT_OUT_LOOSE -o ${X_DOMINANT_OUT_LOOSE}_annotated.txt

# X De Novo
gemini x_linked_de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	-d $LOOSE_MIN_DP \
	--min-gq $LOOSE_MIN_GQ \
	$GEMINIDB > $X_DENOVO_OUT_LOOSE
python $TableAnnotator -i $X_DENOVO_OUT_LOOSE -o ${X_DENOVO_OUT_LOOSE}_annotated.txt

# X Recessive
gemini x_linked_recessive \
        --columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"\
        -d $LOOSE_MIN_DP \
        --min-gq $LOOSE_MIN_GQ \
        $GEMINIDB > $X_RECESSIVE_OUT_LOOSE
python $TableAnnotator -i $X_RECESSIVE_OUT_LOOSE -o ${X_RECESSIVE_OUT_LOOSE}_annotated.txt

# Autosomal Dominant
gemini autosomal_dominant \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"\
	-d $LOOSE_MIN_DP \
	--min-gq $LOOSE_MIN_GQ \
	$GEMINIDB > $AUTODOM_OUT_LOOSE
python $TableAnnotator -i $AUTODOM_OUT_LOOSE -o ${AUTODOM_OUT_LOOSE}_annotated.txt

##########################################################################################

#####################
#       Phase 3     #
#####################

# General Queries to the database:
GNOMAD_GENOME_COMMON='aaf_gnomad_genome_all > 0.01'
GENEHANCER="genehancer is not NULL"
DBSNP="rs_ids is not NULL"


# Check your sample info
gemini query -q "SELECT * FROM samples" $GEMINIDB > $DBSTATFILE


#I will print out a table
echo "Annotation	Total	Common	Rare" >> $DBSTATFILE

# Total number of variants in the file
TOTAL_VARS=`gemini query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP" $GEMINIDB`
rareTOTAL_VARS=`gemini query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP AND $GNOMAD_RARE" $GEMINIDB`
commonTOTAL_VARS=`gemini query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP AND $GNOMAD_COMMON" $GEMINIDB`

echo "TOTAL	$TOTAL_VARS	$commonTOTAL_VARS	$rareTOTAL_VARS" >> $DBSTATFILE

# Exonic variants
EXONIC_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $SEGDUP" $GEMINIDB`
rareEXONIC_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
commonEXONIC_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`

echo "EXONIC	$EXONIC_VARS	$commonEXONIC_VARS	$rareEXONIC_VARS" >> $DBSTATFILE

# Coding Variants
CODING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $SEGDUP" $GEMINIDB`
rareCODING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
commonCODING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`

echo "CODING	$CODING_VARS	$commonCODING_VARS	$rareCODING_VARS" >> $DBSTATFILE

#Splicing variants 
SPLICING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $SEGDUP" $GEMINIDB`
rareSPLICING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
commonSPLICING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`

echo "SPLICING	$SPLICING_VARS	$commonSPLICING_VARS	$rareSPLICING_VARS" >> $DBSTATFILE

# LOF variants
LOF_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $SEGDUP" $GEMINIDB`
rareLOF_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
commonLOF_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`

echo "LOF	$LOF_VARS	$commonLOF_VARS	$rareLOF_VARS" >> $DBSTATFILE

# GeneHancer Variants
GENEHANCER_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $GENEHANCER AND $SEGDUP" $GEMINIDB`
rareGENEHANCER_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $GENEHANCER AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
commonGENEHANCER_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
	WHERE $GENEHANCER AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`

echo "GENEHANCER	$GENEHANCER_VARS	$commonGENEHANCER_VARS	$rareGENEHANCER_VARS" >> $DBSTATFILE


# NOTE: Not much use for this, just takes up space at the moment
## Repeat Above, but use a genotype quality cutoff
#echo "Filtered:" >> $DBSTATFILE
##I will print out a table
#echo "Annotation	Total	Common	Rare" >> $DBSTATFILE
#
## Total number of variants in the file
#TOTAL_VARS=`gemini query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP" $GEMINIDB`
#rareTOTAL_VARS=`gemini query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP AND $GNOMAD_RARE" $GEMINIDB`
#commonTOTAL_VARS=`gemini query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP AND $GNOMAD_COMMON" $GEMINIDB`
#
#echo "TOTAL	$TOTAL_VARS	$commonTOTAL_VARS	$rareTOTAL_VARS" >> $DBSTATFILE
#
## Exonic variants
#EXONIC_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $EXONIC AND $SEGDUP" $GEMINIDB`
#rareEXONIC_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $EXONIC AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
#commonEXONIC_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $EXONIC AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`
#
#echo "EXONIC	$EXONIC_VARS	$commonEXONIC_VARS	$rareEXONIC_VARS" >> $DBSTATFILE
#
## Coding Variants
#CODING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $CODING AND $SEGDUP" $GEMINIDB`
#rareCODING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $CODING AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
#commonCODING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $CODING AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`
#
#echo "CODING	$CODING_VARS	$commonCODING_VARS	$rareCODING_VARS" >> $DBSTATFILE
#
##Splicing variants 
#SPLICING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $SPLICING AND $SEGDUP" $GEMINIDB`
#rareSPLICING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $SPLICING AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
#commonSPLICING_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $SPLICING AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`
#
#echo "SPLICING	$SPLICING_VARS	$commonSPLICING_VARS	$rareSPLICING_VARS" >> $DBSTATFILE
#
## LOF variants
#LOF_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $LOF AND $SEGDUP" $GEMINIDB`
#rareLOF_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $LOF AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
#commonLOF_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $LOF AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`
#
#echo "LOF	$LOF_VARS	$commonLOF_VARS	$rareLOF_VARS" >> $DBSTATFILE
#
## GeneHancer Variants
#GENEHANCER_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $GENEHANCER AND $SEGDUP" $GEMINIDB`
#rareGENEHANCER_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $GENEHANCER AND $SEGDUP AND $GNOMAD_RARE"  $GEMINIDB`
#commonGENEHANCER_VARS=`gemini query -q "SELECT COUNT(*) FROM variants \
#	WHERE $GENEHANCER AND $SEGDUP AND $GNOMAD_COMMON"  $GEMINIDB`
#
#echo "GENEHANCER	$GENEHANCER_VARS	$commonGENEHANCER_VARS	$rareGENEHANCER_VARS" >> $DBSTATFILE
#



#####################
#       Phase 4     #
#####################

# General Damaging Variant List for the proband, which will include clinvar pathogenic hits no matter what
gemini query -q "SELECT $COLUMNS FROM variants WHERE (clinvar_pathogenic is not NULL) AND (clinvar_disease_name is not NULL)" \
	--gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all)" \
	--header \
	$GEMINIDB > $CLINVAR_HITS

python  $TableAnnotator -i $CLINVAR_HITS -o ${CLINVAR_HITS}_annotated.txt


gemini query -q "SELECT $COLUMNS FROM variants WHERE $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	--gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all)" \
	--header \
	$GEMINIDB > $GENERAL_DAMAGING

python $TableAnnotator -i $GENERAL_DAMAGING -o ${GENERAL_DAMAGING}_annotated.txt

gemini query -q "SELECT $COLUMNS FROM variants WHERE $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)" \
	--gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all)" \
	--header \
	$GEMINIDB > $GENERAL_DAMAGING_LOOSE

python $TableAnnotator -i $GENERAL_DAMAGING_LOOSE -o ${GENERAL_DAMAGING_LOOSE}_annotated.txt



#####################
#       Phase 5     #
#####################



# Integrate into single file
# Header 
CVL=$WORKING_DIR${FAMILY_ID}_CVL.tsv
DATE=`date`
echo "FAMILY:	$FAMILY_ID" > $CVL
echo "DATE CVL PRODUCED:	$DATE" >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL

# General Stats
cat $DBSTATFILE >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL

# De novo
echo "De Novo (Strict)" >> $CVL 
cat ${DENOVO_OUT}_annotated.txt >> $CVL


echo "De Novo (Loose)" >> $CVL
cat ${DENOVO_OUT_LOOSE}_annotated.txt >> $CVL

echo "---------------------------------------------------------------------------------------------------------------------" >> $CVL

echo "X De Novo (Strict)" >> $CVL
cat ${X_DENOVO_OUT}_annotated.txt >> $CVL

echo "X De Novo (Loose)" >> $CVL
cat ${X_DENOVO_OUT_LOOSE}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL

# Compound Het
echo "Compound Het (Strict)" >> $CVL
cat ${COMPOUND_HET_OUT}_annotated.txt >> $CVL

echo "Compound Het (Loose)" >> $CVL
cat ${COMPOUND_HET_OUT_LOOSE}_annotated.txt >> $CVL


echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL

# Recessive
echo "Recessive (Strict)" >> $CVL
cat ${RECESSIVE_OUT}_annotated.txt >> $CVL


echo "Recessive (Loose)" >> $CVL
cat ${RECESSIVE_OUT_LOOSE}_annotated.txt >> $CVL

echo "---------------------------------------------------------------------------------------------------------------------" >> $CVL

echo "X Recessive (Strict)" >> $CVL
cat ${X_RECESSIVE_OUT}_annotated.txt >> $CVL

echo "X Recessive (Loose)" >> $CVL
cat ${X_RECESSIVE_OUT_LOOSE}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL

# Dominant 
echo "Autosomal Dominant (Strict)" >> $CVL
cat ${AUTODOM_OUT}_annotated.txt >> $CVL


echo "Autosomal Dominant (Loose)" >> $CVL
cat ${AUTODOM_OUT_LOOSE}_annotated.txt >> $CVL

echo "---------------------------------------------------------------------------------------------------------------------" >> $CVL

echo "X Autosomal Dominant (Strict)" >> $CVL
cat ${X_DOMINANT_OUT}_annotated.txt >> $CVL

echo "X Autosomal Dominant (Loose)" >> $CVL
cat ${X_DOMINANT_OUT_LOOSE}_annotated.txt >> $CVL


echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL

# General Damaging
echo "General Damaging" >> $CVL
cat ${GENERAL_DAMAGING}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL

echo "General Damaging Loose" >> $CVL
cat ${GENERAL_DAMAGING_LOOSE}_annotated.txt >> $CVL


############

# Clean Up
rm *_annotated.txt
rm $AUTODOM_OUT
rm $DENOVO_OUT
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
rm $GENERAL_DAMAGING
rm $GENERAL_DAMAGING_LOOSE

