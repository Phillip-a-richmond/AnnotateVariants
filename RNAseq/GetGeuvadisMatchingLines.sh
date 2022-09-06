# This script will parse through the directory which has the diversity samples (which also have corresponding 1k Genomes samples),
# and grep for the matching sample within Geuvadis samplesheet.tsv

# Diversity samples (From Mark Bennett June 23 2021):
DIVERSITY_DIR=/mnt/scratch/Public/RICHMOND/EHdn/Diversity_1000G_align/

# Thousand genomes samples (From Mark Bennett February 2021)
GEUVADIS_CSV=/mnt/common/OPEN_DATA/GEUVADIS/igsr_Geuvadis.tsv

# We'll make a new directory to make sample naming easier.
SUBSET_CSV=/mnt/common/OPEN_DATA/GEUVADIS/igsr_Geuvadis_Subset.tsv
# add header
head -n1 $GEUVADIS_CSV > $SUBSET_CSV

# Loop through Diversity Dir, extract sample IDs, grep from larger Geuvadis CSV, append to subset CSV
for filename in $(ls $DIVERSITY_DIR/*json); do
	echo $filename
	IFS='/' read -a array <<< $filename
	SJSON=${array[-1]}
	SAMPLE_ID=${SJSON::-17}
	echo $SAMPLE_ID
	
	# grep for the matching sample line, append to subset.csv
	grep $filename $GEUVADIS_CSV >> $SUBSET_CSV 


done
exit


