# Gemini Query script
# Phillip Richmond
# ChangeLog:
# July 10th 2018, script created with addition of gnomad exome and revised filters
# August 30th 2018, added PrimateAI to query, fixed ClinVar to be split to pathogenic/likely pathogenic, fixed General Damaging to be split with Het/Homo
# September 10th, 2018 - Fixed issue with general damaging which impacted reporting of variants
# NOTES: 
# Updates to this script should accomodate differences in how many individuals are present


## Define variables
WORKING_DIR='/mnt/causes-vnx2/TIDE/PROCESS/CYNTHIA_STRAB/'
GEMINIDB='STRAB.db'
FAMILY_ID='STRAB'
TableAnnotator=/mnt/causes-vnx1/PIPELINES/AnnotateVariants/TableAnnotators/GeminiTable2CVL.py

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
GENERAL_DAMAGING_HET=${FAMILY_ID}_GeneralDamaging_Het
GENERAL_DAMAGING_HOMO=${FAMILY_ID}_GeneralDamaging_Homo
CLINVAR_HITS=${FAMILY_ID}_Clinvar_Hits

# Columns to present within CVL
COLUMNS="chrom, start, end, ref, alt, gene, exon, aa_change, impact, impact_severity, rs_ids, filter, gts, gt_ref_depths, gt_alt_depths, gt_quals, in_segdup, confidentregion, inhousedb_ac, gnomad_exome_af_global, gnomad_exome_hom_global, gnomad_exome_af_afr, gnomad_exome_hom_afr, gnomad_exome_af_amr, gnomad_exome_hom_amr, gnomad_exome_af_asj, gnomad_exome_hom_asj, gnomad_exome_af_eas, gnomad_exome_hom_eas, gnomad_exome_af_fin, gnomad_exome_hom_fin, gnomad_exome_af_nfe, gnomad_exome_hom_nfe, gnomad_exome_af_oth, gnomad_exome_hom_oth, gnomad_genome_af_global, gnomad_genome_hom_global, gnomad_genome_af_afr, gnomad_genome_hom_afr, gnomad_genome_af_amr, gnomad_genome_hom_amr, gnomad_genome_af_asj, gnomad_genome_hom_asj, gnomad_genome_af_eas, gnomad_genome_hom_eas, gnomad_genome_af_fin, gnomad_genome_hom_fin, gnomad_genome_af_nfe, gnomad_genome_hom_nfe, gnomad_genome_af_oth, gnomad_genome_hom_oth, primateai, cadd, pp2hdiv, pp2hvar, clinvar_disease_name, clinvar_pathogenic, clinvar_dbinfo"

# Define Variant Annotation Cutoffs
GNOMAD_GENOME_RARE="( (gnomad_genome_af_global <= 0.01) or (gnomad_genome_af_global is NULL)) AND ((gnomad_genome_hom_global <= 10) or (gnomad_genome_hom_global is NULL) )"
GNOMAD_EXOME_RARE="( (gnomad_exome_af_global <= 0.01 ) or (gnomad_exome_af_global is NULL)) AND ( (gnomad_exome_hom_global <= 10) or (gnomad_exome_hom_global is NULL) )"
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

STRICT_MIN_DP=20
STRICT_MIN_GQ=30

LOOSE_MIN_DP=15
LOOSE_MIN_GQ=20

# STRICT FILTER
# --filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $CONFIDENTREGION AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED) AND $FILTER"\

# LOOSE FILTER
# --filter "$GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $INHOUSE_RARE AND ($IMPACT_HIGH OR $IMPACT_MED)"

#########################################################################################

# All variants
gemini query --region 14:20000000-35000000 \
        -q "select $COLUMNS from variants where $CONFIDENTREGION AND $SEGDUP AND $FILTER" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) and (gt_types).(phenotype==1).(==HOM_REF).(all)" \
        --header \
        $GEMINIDB \
        > Strab_Region_Total_All.tsv

python $TableAnnotator -i Strab_Region_Total_All.tsv -o Strab_Region_Total_All.annotated.tsv


# All variants, rare
gemini query --region 14:20000000-35000000 \
	-q "select $COLUMNS from variants where $CONFIDENTREGION AND $SEGDUP AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $FILTER" \
	--gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) and (gt_types).(phenotype==1).(==HOM_REF).(all)" \
	--header \
	$GEMINIDB \
	> Strab_Region_Total_Rare.tsv

python $TableAnnotator -i Strab_Region_Total_Rare.tsv -o Strab_Region_Total_Rare.annotated.tsv


# HIGH / MED Impact 
gemini query --region 14:20000000-35000000 \
        -q "select $COLUMNS from variants where $CONFIDENTREGION AND $SEGDUP AND $FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) and (gt_types).(phenotype==1).(==HOM_REF).(all)" \
        --header \
        $GEMINIDB \
        > Strab_Region_HIGHMED_All.tsv

python $TableAnnotator -i Strab_Region_HIGHMED_All.tsv -o Strab_Region_HIGHMED_All.annotated.tsv

# HIGH / MED Impact Rare
gemini query --region 14:20000000-35000000 \
        -q "select $COLUMNS from variants where $CONFIDENTREGION AND $SEGDUP AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $FILTER AND ($IMPACT_HIGH OR $IMPACT_MED)" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) and (gt_types).(phenotype==1).(==HOM_REF).(all)" \
        --header \
        $GEMINIDB \
        > Strab_Region_HIGHMED_Rare.tsv

python $TableAnnotator -i Strab_Region_HIGHMED_Rare.tsv -o Strab_Region_HIGHMED_Rare.annotated.tsv


# LOW impact
gemini query --region 14:20000000-35000000 \
        -q "select $COLUMNS from variants where $CONFIDENTREGION AND $SEGDUP AND $FILTER AND $IMPACT_LOW" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) and (gt_types).(phenotype==1).(==HOM_REF).(all)" \
        --header \
        $GEMINIDB \
        > Strab_Region_LOW_All.tsv

python $TableAnnotator -i Strab_Region_LOW_All.tsv -o Strab_Region_LOW_All.annotated.tsv

# LOW impact
gemini query --region 14:20000000-35000000 \
        -q "select $COLUMNS from variants where $CONFIDENTREGION AND $SEGDUP AND $GNOMAD_GENOME_RARE AND $GNOMAD_EXOME_RARE AND $FILTER AND $IMPACT_LOW" \
        --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) and (gt_types).(phenotype==1).(==HOM_REF).(all)" \
        --header \
        $GEMINIDB \
        > Strab_Region_LOW_Rare.tsv

python $TableAnnotator -i Strab_Region_LOW_Rare.tsv -o Strab_Region_LOW_Rare.annotated.tsv


