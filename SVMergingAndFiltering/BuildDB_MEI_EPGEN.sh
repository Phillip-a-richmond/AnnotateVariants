# Build MEI EPGEN
# Author: Phillip Richmond
# Contact: phillip.a.richmond@gmail.com
# This pipeline is used to create a merged database (BED format) for a set of MEI VCFs which were called with MELT.
# This specific analysis is for the EPGEN data, which has mergedMEI.vcf files from joint sample calling with MELT.
# The flow here is: 
# 0-setup
# 1-split joined VCFs into individual sample vcfs
# 2-merge the individual sample VCFs with SURVIVOR
# 3-convert from VCF to bedpe with SURVIVOR
# 4-add frequency from VCF to bedpe with a little python script
# 5-cut from bedpe to bed
# 6-add SV name to bed
# 7-cutout count and frequency into mini bed files and header them for use with AnnotSV



# Step 0 - Setup
## Define BuildInHouseDB dir
BUILDINHOUSEDBDIR=/mnt/common/Precision/BuildInHouseDB/
## Load executables
source $BUILDINHOUSEDBDIR/opt/miniconda3/etc/profile.d/conda.sh
conda activate $BUILDINHOUSEDBDIR/opt/InHouseDB_environment/

## SURVIVOR executable, if conda works then this is in your path now
SURVIVOR=SURVIVOR

## BuildInhouseDB dir
BUILDINHOUSEDB=/mnt/common/Precision/BuildInHouseDB/

## Set the input VCF directory
VCFDIR=/mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/

## Input data files
MEIFILES=$VCFDIR/*MergedMEI.vcf

## AnnotSV dir
ANNOTSVDIR=/mnt/common/Precision/AnnotSV/share/AnnotSV/Annotations_Human/Users/GRCh38/



# Step 1 - Split each VCF to it's subfile components
## Split Family MERGED VCF to get sub-VCFs for each individual

## Split MEI files
for vcf in $MEIFILES
do
	echo $vcf
	bgzip -c $vcf > ${vcf}.gz
	tabix ${vcf}.gz
	file=${vcf}.gz
	for sample in `bcftools query -l $file`; do
	    bcftools view -c1 -Ov -s $sample -o ${file/.vcf*/.$sample.vcf} $file
	done
	
done


# Step 2 - merge the VCFs

## Get the date for tracking build purposes
DateOfCreation=`date +%Y%m%d`

## Get the list of VCFs into a simple file
## Here we're using the specific VCFs from EPGEN, but change this for other sample sets
cd $VCFDIR
ls *MergedMEI.EPGEN*.vcf > Sample_MEI_VCFs_list_${DateOfCreation}.txt

## SURVIVOR merge command
# File with VCF names and paths
# max distance between breakpoints (0-1 percent of length, 1- number of bp) 
# Minimum number of supporting caller
# Take the type into account (1==yes, else no)
# Take the strands of SVs into account (1==yes, else no)
# Disabled.
# Minimum size of SVs to be taken into account.
# Output VCF filename
$SURVIVOR merge \
	 Sample_MEI_VCFs_list_${DateOfCreation}.txt \
	 25 \
	 1 \
	 1 \
	 0 \
	 0 \
	 1 \
	 InHouseDB_MEI_15_${DateOfCreation}.vcf


# Step 3 - From VCF --> BEDPE
## SURVIVOR vcftobed
# vcf file
# min size
# max size
# output file
$SURVIVOR vcftobed \
	InHouseDB_MEI_15_${DateOfCreation}.vcf \
	0 \
	1000000000 \
	InHouseDB_MEI_15_${DateOfCreation}.bedpe



# Step 4 - Add frequency to bedpe
python $BUILDINHOUSEDBDIR/AddFreqToBedpe.py -V InHouseDB_MEI_15_${DateOfCreation}.vcf \
	-B InHouseDB_MEI_15_${DateOfCreation}.bedpe \
	-O InHouseDB_MEI_15_${DateOfCreation}.freq.bedpe


# Step 5 - Cut out columns and make a bed file
# Also remove TRA
grep -v 'TRA' InHouseDB_MEI_15_${DateOfCreation}.freq.bedpe | \
	cut -f1,2,6,7,11,12,13 - \
	> InHouseDB_MEI_15_${DateOfCreation}.bed

# Step 6 - add SV name 
python $BUILDINHOUSEDBDIR/add_SV_name.py \
	InHouseDB_MEI_15_${DateOfCreation}.bed \
	InHouseDB_MEI_15_${DateOfCreation}.final.bed


# Step 7 - add header line, cut-out freq and af columns from final bed, concatenate with header
## Frequency 
echo "#chrom	start	end	InHouseDB_MEI_Freq" > InHouseDB_MEI_Freq.header
cut -f1,2,3,7  InHouseDB_MEI_15_${DateOfCreation}.final.bed | cat InHouseDB_MEI_Freq.header - > UserInHouseDB_MEI_${DateOfCreation}.Freq.bed


## Count
echo "#chrom	start	end	InHouseDB_MEI_Count" > InHouseDB_MEI_Count.header
cut -f1,2,3,6  InHouseDB_MEI_15_${DateOfCreation}.final.bed | cat InHouseDB_MEI_Count.header - > UserInHouseDB_MEI_${DateOfCreation}.Count.bed

# Step 8 - Copy files to AnnotSV folder
cp UserInHouseDB* $ANNOTSVDIR/SVincludedInFt/
cp UserInHouseDB* $ANNOTSVDIR/FtIncludedInSV/
cp UserInHouseDB* $ANNOTSVDIR/AnyOverlap/
exit


#######
# END #
#######






## Old code that the above was motivated by. 

## Splitting at the end could still be useful, or iterating over different overlap distance threshholds


# Running for both max_dist between edges of 100, and 1000
for MaxDist in {50,100}
do

	# These commands to SURVIVOR invoke this:
	# mergevcf: Consensus Call from multiple SV vcf files
	# Sample_VCFs_list_${DateOfCreation}.txt - Tab file with names
	# max distance between breakpoints - The $MaxDist variable
	# Minimum number of supporting caller - 1
	# Take the type into account (1==yes, else no) - 1
	#  Take the strands of SVs into account (1==yes, else no) - 0
	# Estimate distance based on the size of SV (1==yes, else no). -1
	# Minimum size of SVs to be taken into account. - 50
	# Output prefix 
	
	# Run SURVIVOR
	$SURVIVOR merge \
		Sample_VCFs_list_${DateOfCreation}.txt \
		$MaxDist \
		1 \
		1 \
		0 \
	 	1 \
		50 \
		InHouseDB_SV_${MaxDist}_${DateOfCreation}.vcf
	
	# For Visualization purposes, sort and index the vcf file.
	igvtools sort InHouseDB_SV_${MaxDist}_${DateOfCreation}.vcf  InHouseDB_SV_${MaxDist}_${DateOfCreation}_sorted.vcf
	igvtools index InHouseDB_SV_${MaxDist}_${DateOfCreation}_sorted.vcf
	
	
	# These commands to SURVIVOR invoke this with option 8:
	# vcftobedpe: Convert vcf to bedpe
	# vcf file -InHouseDB_SV_${MaxDist}_${DateOfCreation}.vcf
	# min size
	# max size - 15000000, I set this to 15MB just to get rid of those stupid like half-chromosome artifacts or fully deleted centromere artifacts, Could easily increase this
	# output file - InHouseDB_SV_${MaxDist}_${DateOfCreation}.bedpe
	
	# Convert into BEDPE file 
	$SURVIVOR vcftobed \
		InHouseDB_SV_${MaxDist}_${DateOfCreation}.vcf \
		1 \
		5000000 \
		InHouseDB_SV_${MaxDist}_${DateOfCreation}.bedpe	
	
	
	# Use the BEDPE and the VCF in order to add the SUPP=# to the BEDPE
	python AddFreqToBedpe.py -V InHouseDB_SV_${MaxDist}_${DateOfCreation}.vcf -B InHouseDB_SV_${MaxDist}_${DateOfCreation}.bedpe -O InHouseDB_SV_${MaxDist}_${DateOfCreation}_counted.bedpe
	
	# Cut out columns and make a bed file
	cut -f1,2,6,7,11,12,13 InHouseDB_SV_${MaxDist}_${DateOfCreation}_counted.bedpe > InHouseDB_SV_${MaxDist}_${DateOfCreation}.bed
	
	grep "DEL" InHouseDB_SV_${MaxDist}_${DateOfCreation}.bed > InHouseDB_SV_${MaxDist}_${DateOfCreation}_DEL.bed
	grep "DUP" InHouseDB_SV_${MaxDist}_${DateOfCreation}.bed > InHouseDB_SV_${MaxDist}_${DateOfCreation}_DUP.bed
	grep "INS" InHouseDB_SV_${MaxDist}_${DateOfCreation}.bed > InHouseDB_SV_${MaxDist}_${DateOfCreation}_INS.bed
	grep "INV" InHouseDB_SV_${MaxDist}_${DateOfCreation}.bed > InHouseDB_SV_${MaxDist}_${DateOfCreation}_INV.bed

	python add_SV_name.py InHouseDB_SV_${MaxDist}_${DateOfCreation}_DEL.bed InHouseDB_SV_${MaxDist}_${DateOfCreation}_DEL_name.bed
	python add_SV_name.py InHouseDB_SV_${MaxDist}_${DateOfCreation}_DUP.bed InHouseDB_SV_${MaxDist}_${DateOfCreation}_DUP_name.bed
	python add_SV_name.py InHouseDB_SV_${MaxDist}_${DateOfCreation}_INV.bed InHouseDB_SV_${MaxDist}_${DateOfCreation}_INV_name.bed
	python add_SV_name.py InHouseDB_SV_${MaxDist}_${DateOfCreation}_INS.bed InHouseDB_SV_${MaxDist}_${DateOfCreation}_INS_name.bed
		
done	










