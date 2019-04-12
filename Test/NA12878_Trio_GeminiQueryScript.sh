WORKING_DIR=/mnt/causes-vnx1/PIPELINES/AnnotateVariants/Test/
GEMINIDB=NA12878_Trio.db
FAMILY_ID=NA12878_Trio
TableAnnotator=/mnt/causes-vnx1/PIPELINES/AnnotateVariants//TableAnnotators/GeminiTable2CVL-devReOrder.py

# ^^ Header above describes family, location of GeminiDB, location of table annotator 

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
	

source /opt/tools/hpcenv.sh
GEMINI=/mnt/causes-vnx1/DATABASES/GEMINI-2019/bin/gemini
#GEMINI=gemini
PATH=/mnt/causes-vnx1/DATABASES/GEMINI-2019/:$PATH
## Change to working directory
cd $WORKING_DIR

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
AUTODOM_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_autoDom_loose
DENOVO_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_deNovo_loose
#DENOVO_LOW_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_deNovoLow_loose; NOTE: This is too long of a list, removing for now
RECESSIVE_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_recessive_loose
COMPOUND_HET_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_compoundHet_loose
X_RECESSIVE_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_Xrecessive_loose
X_DOMINANT_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_Xdominant_loose
X_DENOVO_OUT_LOOSE=$WORKING_DIR${FAMILY_ID}_Xdenovo_loose

# Phase 3
DBSTATFILE=${FAMILY_ID}_GeminiDB_Stats.txt

# Phase 4
GENERAL_DAMAGING_HET=${FAMILY_ID}_GeneralDamaging_Het
GENERAL_DAMAGING_HOMO=${FAMILY_ID}_GeneralDamaging_Homo
CLINVAR_HITS=${FAMILY_ID}_Clinvar_Hits
SPLICING_HITS=${FAMILY_ID}_SpliceCandidates
NONCODING_HITS=${FAMILY_ID}_NoncodingCandidates

# Columns to present within CVL
# OLD COLUMNS="chrom, start, end, ref, alt, gene, exon, aa_change, impact, impact_severity, rs_ids, filter, gts, gt_ref_depths, gt_alt_depths, gt_quals, in_segdup, confidentregion, inhousedb_ac, gnomad_exome_af_global, gnomad_exome_hom_global, gnomad_exome_af_afr, gnomad_exome_hom_afr, gnomad_exome_af_amr, gnomad_exome_hom_amr, gnomad_exome_af_asj, gnomad_exome_hom_asj, gnomad_exome_af_eas, gnomad_exome_hom_eas, gnomad_exome_af_fin, gnomad_exome_hom_fin, gnomad_exome_af_nfe, gnomad_exome_hom_nfe, gnomad_exome_af_oth, gnomad_exome_hom_oth, gnomad_genome_af_global, gnomad_genome_hom_global, gnomad_genome_af_afr, gnomad_genome_hom_afr, gnomad_genome_af_amr, gnomad_genome_hom_amr, gnomad_genome_af_asj, gnomad_genome_hom_asj, gnomad_genome_af_eas, gnomad_genome_hom_eas, gnomad_genome_af_fin, gnomad_genome_hom_fin, gnomad_genome_af_nfe, gnomad_genome_hom_nfe, gnomad_genome_af_oth, gnomad_genome_hom_oth, primateai, cadd, pp2hdiv, pp2hvar, clinvar_disease_name, clinvar_pathogenic, clinvar_dbinfo"
COLUMNS="chrom, start, end, ref, alt, gene, exon, aa_change, codon_change, transcript, biotype, impact, impact_severity, rs_ids, filter, gts, gt_ref_depths, gt_alt_depths, gt_alt_freqs, gt_quals, in_segdup, confidentregion, inhousedb_ac, gnomad_exome_ac_global, gnomad_exome_af_global, gnomad_exome_hom_global, gnomad_exome_popmax, gnomad_exome_af_popmax, gnomad_exome_hom_popmax, gnomad_exome_AF_controls, gnomad_exome_hom_controls, gnomad_genome_ac_global, gnomad_genome_af_global, gnomad_genome_hom_global, gnomad_genome_popmax, gnomad_genome_af_popmax, gnomad_genome_hom_popmax, gnomad_genome_AF_controls, gnomad_genome_hom_controls, primateai, cadd, cadd_indel, pp2hdiv, pp2hvar, spliceai_acceptorgain, spliceai_acceptorloss, spliceai_donorgain, spliceai_donorloss, fathmm_xf_noncoding, fitcons, genehancer, clinvar_disease_name, clinvar_pathogenic, clinvar_dbinfo, ccr"

# Define Variant Annotation Cutoffs
GNOMAD_GENOME_RARE="( (gnomad_genome_af_global <= 0.001) or (gnomad_genome_af_global is NULL)) AND ((gnomad_genome_hom_global <= 10) or (gnomad_genome_hom_global is NULL) )"
GNOMAD_EXOME_RARE="( (gnomad_exome_af_global <= 0.001 ) or (gnomad_exome_af_global is NULL)) AND ( (gnomad_exome_hom_global <= 10) or (gnomad_exome_hom_global is NULL) )"
GNOMAD_GENOME_DENOVO_RARE="( (gnomad_genome_ac_global <= 5) or (gnomad_genome_af_global is NULL)) AND ((gnomad_genome_hom_global <= 10) or (gnomad_genome_hom_global is NULL) )"
GNOMAD_EXOME_DENOVO_RARE="( (gnomad_exome_ac_global <= 5) or (gnomad_exome_af_global is NULL)) AND ( (gnomad_exome_hom_global <= 10) or (gnomad_exome_hom_global is NULL) )"
INHOUSE_RARE='(inhousedb_ac <= 3 or inhousedb_ac is NULL)'

CODING='is_coding=1'
EXONIC='is_exonic=1'
SPLICING='is_splicing=1'

LOF='is_lof=1'

IMPACT_HIGH="impact_severity=='HIGH'"
IMPACT_MED="impact_severity=='MED'"
IMPACT_LOW="impact_severity=='LOW'"

FILTER='filter is NULL'
UNFILTER='filter is not NULL'
CONFIDENTREGION='confidentregion = 1'
SEGDUP='in_segdup=0'

# Chose 0.5 based on spliceAI manuscript
# https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf
SPLICEAI='((spliceai_acceptorgain >= 0.5) OR (spliceai_acceptorloss >= 0.5 ) OR (spliceai_donorgain >= 0.5) OR (spliceai_donorloss >= 0.5 ))'

CADD='((cadd >= 20) OR (cadd_indel >= 20))'
FATHMM_NONCODING='(fathmm_xf_noncoding >= 0.9)'
NONCODING="NOT $CODING AND ($FATHMM_NONCODING OR $CADD)"

STRICT_MIN_DP=15
STRICT_MIN_GQ=30

# Dropping in pipeline update 20190321
#LOOSE_MIN_DP=15
#LOOSE_MIN_GQ=20

# STRICT FILTER
# --filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\

# LOOSE FILTER
# --filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"

#########################################################################################

#####################
#       Phase 1     #
#####################
# Inheritance models
## Recessive
$GEMINI autosomal_recessive \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $RECESSIVE_OUT
python $TableAnnotator -i $RECESSIVE_OUT -o ${RECESSIVE_OUT}_annotated.txt

# Compound Het Variants
$GEMINI comp_hets \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $COMPOUND_HET_OUT
python $TableAnnotator -i $COMPOUND_HET_OUT -o ${COMPOUND_HET_OUT}_annotated.txt

# De novo
$GEMINI de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $DENOVO_OUT
python $TableAnnotator -i $DENOVO_OUT -o ${DENOVO_OUT}_annotated.txt

# De novo Not HIGH/MED
$GEMINI de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_DENOVO_RARE AND $GNOMAD_EXOME_DENOVO_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND NOT ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $DENOVO_LOW_OUT
python $TableAnnotator -i $DENOVO_LOW_OUT -o ${DENOVO_LOW_OUT}_annotated.txt

# X Dominant
$GEMINI x_linked_dominant  \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $X_DOMINANT_OUT
python $TableAnnotator -i $X_DOMINANT_OUT -o ${X_DOMINANT_OUT}_annotated.txt

# X De Novo
$GEMINI x_linked_de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $X_DENOVO_OUT
python $TableAnnotator -i $X_DENOVO_OUT -o ${X_DENOVO_OUT}_annotated.txt

# X Recessive
$GEMINI x_linked_recessive \
        --columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
        -d $STRICT_MIN_DP \
        --min-gq $STRICT_MIN_GQ \
        $GEMINIDB > $X_RECESSIVE_OUT
python $TableAnnotator -i $X_RECESSIVE_OUT -o ${X_RECESSIVE_OUT}_annotated.txt

# Autosomal Dominant
$GEMINI autosomal_dominant \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $AUTODOM_OUT
python $TableAnnotator -i $AUTODOM_OUT -o ${AUTODOM_OUT}_annotated.txt

# Mendel Errors
$GEMINI mendel_errors \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	-d $STRICT_MIN_DP \
	--min-gq $STRICT_MIN_GQ \
	$GEMINIDB > $MENDEL_ERRORS_OUT
python $TableAnnotator -i $MENDEL_ERRORS_OUT -o ${MENDEL_ERRORS_OUT}_annotated.txt

##########################################################################################

#####################
#       Phase 2     #
#####################
# Repeat but with Looser filters, specifically on genotype quality, depth, and whether or not it's in a segdup

# Inheritance models
## Recessive
$GEMINI autosomal_recessive \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)  AND NOT ($CONFIDENTREGION AND $SEGDUP AND $FILTER) "\
	$GEMINIDB > $RECESSIVE_OUT_LOOSE
python $TableAnnotator -i $RECESSIVE_OUT_LOOSE -o ${RECESSIVE_OUT_LOOSE}_annotated.txt


# Compound Het Variants
$GEMINI comp_hets \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND NOT ($CONFIDENTREGION AND $SEGDUP AND $FILTER) "\
	$GEMINIDB > $COMPOUND_HET_OUT_LOOSE
python $TableAnnotator -i $COMPOUND_HET_OUT_LOOSE -o ${COMPOUND_HET_OUT_LOOSE}_annotated.txt

# De novo
$GEMINI de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND NOT ($CONFIDENTREGION AND $SEGDUP AND $FILTER) "\
	$GEMINIDB > $DENOVO_OUT_LOOSE
python $TableAnnotator -i $DENOVO_OUT_LOOSE -o ${DENOVO_OUT_LOOSE}_annotated.txt

# De novo LOW
$GEMINI de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND NOT ($IMPACT_HIGH OR $IMPACT_MED) AND NOT ($CONFIDENTREGION AND $SEGDUP AND $FILTER) "\
	$GEMINIDB > $DENOVO_LOW_OUT_LOOSE
python $TableAnnotator -i $DENOVO_OUT_LOOSE -o ${DENOVO_LOW_OUT_LOOSE}_annotated.txt

# X Dominant
$GEMINI x_linked_dominant  \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND NOT ($CONFIDENTREGION AND $SEGDUP AND $FILTER)"\
	$GEMINIDB > $X_DOMINANT_OUT_LOOSE
python $TableAnnotator -i $X_DOMINANT_OUT_LOOSE -o ${X_DOMINANT_OUT_LOOSE}_annotated.txt

# X De Novo
$GEMINI x_linked_de_novo \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND NOT ($CONFIDENTREGION AND $SEGDUP AND $FILTER) "\
	$GEMINIDB > $X_DENOVO_OUT_LOOSE
python $TableAnnotator -i $X_DENOVO_OUT_LOOSE -o ${X_DENOVO_OUT_LOOSE}_annotated.txt

# X Recessive
$GEMINI x_linked_recessive \
        --columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND NOT ($CONFIDENTREGION AND $SEGDUP AND $FILTER) "\
        $GEMINIDB > $X_RECESSIVE_OUT_LOOSE
python $TableAnnotator -i $X_RECESSIVE_OUT_LOOSE -o ${X_RECESSIVE_OUT_LOOSE}_annotated.txt

# Autosomal Dominant
$GEMINI autosomal_dominant \
	--columns "$COLUMNS" \
	--filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND NOT ($CONFIDENTREGION AND $SEGDUP AND $FILTER) "\
	$GEMINIDB > $AUTODOM_OUT_LOOSE
python $TableAnnotator -i $AUTODOM_OUT_LOOSE -o ${AUTODOM_OUT_LOOSE}_annotated.txt

##########################################################################################

#####################
#       Phase 3     #
#####################

# General Queries to the database:
GNOMAD_GENOME_COMMON='(gnomad_genome_af_global > 0.001)'
GNOMAD_EXOME_COMMON='(gnomad_exome_af_global > 0.001)'
GNOMAD_GENOME_RARE='((gnomad_genome_af_global <= 0.001) or (gnomad_genome_af_global is NULL)) '
GNOMAD_EXOME_RARE='((gnomad_exome_af_global <= 0.001) or (gnomad_exome_af_global is NULL))'

GENEHANCER='genehancer is not NULL'
DBSNP='rs_ids is not NULL'
SEGDUP='in_segdup = 0'

UTR='( (is_coding = 0) and (is_exonic = 1))'
# Check your sample info
$GEMINI query -q "SELECT * FROM samples" $GEMINIDB > $DBSTATFILE


#I will print out a table
echo "Annotation	Total	Common	Rare" >> $DBSTATFILE

# Total number of variants in the file
TOTAL_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareTOTAL_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE" $GEMINIDB`
commonTOTAL_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants WHERE $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)" $GEMINIDB`


echo "TOTAL	$TOTAL_VARS	$commonTOTAL_VARS	$rareTOTAL_VARS" >> $DBSTATFILE

# UTR vars
UTR_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $UTR AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareUTR_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $UTR AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonUTR_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $UTR AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "UTR	$UTR_VARS	$commonUTR_VARS	$rareUTR_VARS" >> $DBSTATFILE

# Exonic variants
EXONIC_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareEXONIC_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonEXONIC_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $EXONIC AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "EXONIC	$EXONIC_VARS	$commonEXONIC_VARS	$rareEXONIC_VARS" >> $DBSTATFILE

# Coding Variants
CODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareCODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonCODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $CODING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "CODING	$CODING_VARS	$commonCODING_VARS	$rareCODING_VARS" >> $DBSTATFILE

#Splicing variants 
SPLICING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareSPLICING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonSPLICING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $SPLICING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "SPLICING	$SPLICING_VARS	$commonSPLICING_VARS	$rareSPLICING_VARS" >> $DBSTATFILE

# LOF variants
LOF_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareLOF_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonLOF_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $LOF AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "LOF	$LOF_VARS	$commonLOF_VARS	$rareLOF_VARS" >> $DBSTATFILE

# GeneHancer Variants
GENEHANCER_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $GENEHANCER AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareGENEHANCER_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $GENEHANCER AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonGENEHANCER_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
	WHERE $GENEHANCER AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "GENEHANCER	$GENEHANCER_VARS	$commonGENEHANCER_VARS	$rareGENEHANCER_VARS" >> $DBSTATFILE


# Splicing Variants
SPLICEAI_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $SPLICEAI AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareSPLICEAI_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $SPLICEAI AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonSPLICEAI_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $SPLICEAI AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "SPLICEAI        $SPLICEAI_VARS        $commonSPLICEAI_VARS  $rareSPLICEAI_VARS" >> $DBSTATFILE

# CADD >= 20
CADD_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $CADD AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareCADD_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $CADD AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonCADD_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $CADD AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "CADD        $CADD_VARS        $commonCADD_VARS  $rareCADD_VARS" >> $DBSTATFILE


# FATHMM >= 0.9 (arbitrary, need to figure this out better)
FATHMM_NONCODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $FATHMM_NONCODING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION" $GEMINIDB`
rareFATHMM_NONCODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $FATHMM_NONCODING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE"  $GEMINIDB`
commonFATHMM_NONCODING_VARS=`$GEMINI query -q "SELECT COUNT(*) FROM variants \
        WHERE $FATHMM_NONCODING AND $SEGDUP AND $FILTER AND $CONFIDENTREGION AND ($GNOMAD_GENOME_COMMON OR $GNOMAD_EXOME_COMMON)"  $GEMINIDB`

echo "FATHMM_NONCODING        $FATHMM_NONCODING_VARS        $commonFATHMM_NONCODING_VARS  $rareFATHMM_NONCODING_VARS" >> $DBSTATFILE





#####################
#       Phase 4     #
#####################
#GNOMAD_GENOME_RARE="( (gnomad_genome_af_global <= 0.001 or gnomad_genome_af_global is NULL) AND (gnomad_genome_hom_global <= 10 or gnomad_genome_hom_global is NULL) )"
#GNOMAD_EXOME_RARE="( (gnomad_exome_af_global <= 0.001 or gnomad_exome_af_global is NULL) AND (gnomad_exome_hom_global <= 10 or gnomad_exome_hom_global is NULL) )"
# NOTE: Changing this allows homozygotes to escape into the general damaging. we are fine with that.
GNOMAD_GENOME_RARE='((gnomad_genome_af_global <= 0.001) or (gnomad_genome_af_global is NULL)) '
GNOMAD_EXOME_RARE='((gnomad_exome_af_global <= 0.001) or (gnomad_exome_af_global is NULL))'

# General Damaging Variant List for the proband, which will include clinvar pathogenic hits no matter what
# note: Added rarity to clinvar here
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE ( (clinvar_pathogenic == 'Pathogenic') OR (clinvar_pathogenic == 'Likely Pathogenic') ) AND (clinvar_disease_name is not NULL) AND $CONFIDENTREGION AND $SEGDUP AND $FILTER AND $GNOMAD_EXOME_RARE AND $GNOMAD_GENOME_RARE" \
	--gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(any)" \
	--header \
	$GEMINIDB > $CLINVAR_HITS

python  $TableAnnotator -i $CLINVAR_HITS -o ${CLINVAR_HITS}_annotated.txt

# General damaging het (HIGH/MED impact)
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE $GNOMAD_EXOME_RARE AND $GNOMAD_GENOME_RARE AND $CONFIDENTREGION AND $SEGDUP AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\
	--gt-filter "(gt_types).(phenotype==2).(==HET).(any)" \
	--header \
	$GEMINIDB > $GENERAL_DAMAGING_HET

python $TableAnnotator -i $GENERAL_DAMAGING_HET -o ${GENERAL_DAMAGING_HET}_annotated.txt

# General damaging homo (HIGH/MED impact)
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE $GNOMAD_EXOME_RARE AND $GNOMAD_GENOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $CONFIDENTREGION AND $SEGDUP AND $FILTER" \
	--gt-filter "(gt_types).(phenotype==2).(==HOM_ALT).(any)" \
	--header \
	$GEMINIDB > $GENERAL_DAMAGING_HOMO

python $TableAnnotator -i $GENERAL_DAMAGING_HOMO -o ${GENERAL_DAMAGING_HOMO}_annotated.txt

# Possible Intronic
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE $GNOMAD_EXOME_RARE AND $GNOMAD_GENOME_RARE AND $INHOUSE_RARE AND $SPLICEAI AND $CONFIDENTREGION AND $SEGDUP AND $FILTER" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(any)" \
        --header \
        $GEMINIDB > $SPLICING_HITS
python $TableAnnotator -i $SPLICING_HITS -o ${SPLICING_HITS}_annotated.txt

# Noncoding high impact
$GEMINI query -q "SELECT $COLUMNS FROM variants WHERE $GNOMAD_EXOME_RARE AND $GNOMAD_GENOME_RARE AND $INHOUSE_RARE AND $NONCODING AND $CONFIDENTREGION AND $SEGDUP AND $FILTER" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(any)" \
        --header \
        $GEMINIDB > $NONCODING_HITS
python $TableAnnotator -i $NONCODING_HITS -o ${NONCODING_HITS}_annotated.txt

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
echo "" >> $CVL

# General Stats
cat $DBSTATFILE >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# De novo
echo "De Novo (Strict)" >> $CVL 
cat ${DENOVO_OUT}_annotated.txt >> $CVL

echo "De Novo (Loose)" >> $CVL
cat ${DENOVO_OUT_LOOSE}_annotated.txt >> $CVL

echo "---------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# X de novo
echo "X De Novo (Strict)" >> $CVL
cat ${X_DENOVO_OUT}_annotated.txt >> $CVL

echo "X De Novo (Loose)" >> $CVL
cat ${X_DENOVO_OUT_LOOSE}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Compound Het
echo "Compound Het (Strict)" >> $CVL
cat ${COMPOUND_HET_OUT}_annotated.txt >> $CVL

echo "Compound Het (Loose)" >> $CVL
cat ${COMPOUND_HET_OUT_LOOSE}_annotated.txt >> $CVL


echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Recessive
echo "Recessive (Strict)" >> $CVL
cat ${RECESSIVE_OUT}_annotated.txt >> $CVL


echo "Recessive (Loose)" >> $CVL
cat ${RECESSIVE_OUT_LOOSE}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

echo "X Recessive (Strict)" >> $CVL
cat ${X_RECESSIVE_OUT}_annotated.txt >> $CVL

echo "X Recessive (Loose)" >> $CVL
cat ${X_RECESSIVE_OUT_LOOSE}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Dominant 
echo "Autosomal Dominant (Strict)" >> $CVL
cat ${AUTODOM_OUT}_annotated.txt >> $CVL


echo "Autosomal Dominant (Loose)" >> $CVL
cat ${AUTODOM_OUT_LOOSE}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

echo "X Autosomal Dominant (Strict)" >> $CVL
cat ${X_DOMINANT_OUT}_annotated.txt >> $CVL

echo "X Autosomal Dominant (Loose)" >> $CVL
cat ${X_DOMINANT_OUT_LOOSE}_annotated.txt >> $CVL


echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Clinvar Hits
echo "Clinvar Hits" >> $CVL
cat ${CLINVAR_HITS}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# General Damaging
echo "General Damaging Homozygous" >> $CVL
cat ${GENERAL_DAMAGING_HOMO}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

echo "General Damaging Heterozygous" >> $CVL
cat ${GENERAL_DAMAGING_HET}_annotated.txt >> $CVL

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL


#######################
# Other, deeper dives #
#######################

# De novo low
echo "De Novo Low (Strict)" >> $CVL
cat ${DENOVO_LOW_OUT}_annotated.txt >> $CVL 

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Splicing Candidates 

echo "Splicing Candidates">> $CVL
cat ${SPLICING_HITS}_annotated.txt >> $CVL
echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $CVL
echo "" >> $CVL

# Noncoding Others with CADD / FATHMM-XF
echo "Noncoding Candidates, cis-regulatory"
cat ${NONCODING_HITS}_annotated.txt >> $CVL




############

# Clean Up
#rm *_annotated.txt
#rm $AUTODOM_OUT
#rm $DENOVO_OUT
#rm $DENOVO_LOW_OUT
#rm $RECESSIVE_OUT
#rm $COMPOUND_HET_OUT
#rm $X_RECESSIVE_OUT
#rm $X_DOMINANT_OUT
#rm $X_DENOVO_OUT
#rm $AUTODOM_OUT_LOOSE
#rm $DENOVO_OUT_LOOSE
#rm $RECESSIVE_OUT_LOOSE
#rm $COMPOUND_HET_OUT_LOOSE
#rm $X_RECESSIVE_OUT_LOOSE
#rm $X_DOMINANT_OUT_LOOSE
#rm $X_DENOVO_OUT_LOOSE
#rm $DENOVO_LOW_OUT_LOOSE
#rm $GENERAL_DAMAGING_HET
#rm $GENERAL_DAMAGING_HOMO
#rm $CLINVAR_HITS
#rm $SPLICING_HITS
#rm $NONCODING_HITS
#rm $MENDEL_ERRORS_OUT
