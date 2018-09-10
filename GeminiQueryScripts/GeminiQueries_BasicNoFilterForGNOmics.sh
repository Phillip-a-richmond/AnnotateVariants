## Define variables
WORKING_DIR='/mnt/causes-vnx1/TESTING/'
GEMINIDB='NA12878_Trio.db'
FAMILY_ID='NA12878_Trio'
source /opt/tools/hpcenv.sh

## Change to working directory
cd $WORKING_DIR

# Define Inheritance Model Variant Files
AUTODOM_OUT=$WORKING_DIR${FAMILY_ID}_autoDom
DENOVO_OUT=$WORKING_DIR${FAMILY_ID}_deNovo
RECESSIVE_OUT=$WORKING_DIR${FAMILY_ID}_recessive
COMPOUND_HET_OUT=$WORKING_DIR${FAMILY_ID}_compoundHet
COMPOUND_HET_LOOSE_OUT=$WORKING_DIR${FAMILY_ID}_compoundHet_Loose
X_RECESSIVE_OUT=$WORKING_DIR${FAMILY_ID}_Xrecessive
X_DOMINANT_OUT=$WORKING_DIR${FAMILY_ID}_Xdominant
X_DENOVO_OUT=$WORKING_DIR${FAMILY_ID}_Xdenovo
GENERAL_DAMAGING=${FAMILY_ID}_GeneralDamaging.txt

COLUMNS='variant_id, chrom, start, end, ref, alt, gene'

#####################
#       Phase 1     #
#####################
# Inheritance models
## Recessive
gemini autosomal_recessive \
	--columns "$COLUMNS" \
	$GEMINIDB > $RECESSIVE_OUT

# Compound Het Variants
gemini comp_hets \
	--columns "$COLUMNS" \
	$GEMINIDB > $COMPOUND_HET_OUT

#Compound Het Variants, lower restriction
gemini comp_hets \
	--columns "$COLUMNS" \
	--max-priority 3 \
	$GEMINIDB > $COMPOUND_HET_LOOSE_OUT


# De novo
gemini de_novo \
	--columns "$COLUMNS" \
	$GEMINIDB > $DENOVO_OUT

# X Dominant
gemini x_linked_dominant  \
	--columns "$COLUMNS" \
	$GEMINIDB > $X_DOMINANT_OUT

# X De Novo
gemini x_linked_de_novo \
	--columns "$COLUMNS" \
	$GEMINIDB > $X_DENOVO_OUT

# X Recessive
gemini x_linked_recessive \
        --columns "$COLUMNS" \
        $GEMINIDB > $X_RECESSIVE_OUT

# Autosomal Dominant
gemini autosomal_dominant \
        --columns "$COLUMNS" \
	$GEMINIDB > $AUTODOM_OUT

# General Damaging Variant List for the proband
gemini query -q "SELECT $COLUMNS FROM variants WHERE (impact=='HIGH' OR impact=='MED')" \
	--header \
	$GEMINIDB > $GENERAL_DAMAGING


