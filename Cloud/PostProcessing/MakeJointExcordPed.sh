# This script takes in a directory of excord .bed.gz files, and makes a single file with format:
# sampleID	Filename.Bed.gz

# I'll do that with this code:
shopt -s nullglob
DataDir=/mnt/common/OPEN_DATA/1kG_Trio/Excord_Output/
Files=($DataDir/*.bed.gz)

pedfile=$DataDir/AllSamples.ped
# Delete if it exists
rm $pedfile
# give it a header
echo "Sample	Alt_File" > $pedfile

for file in ${Files[@]};
do
	basefile=$(basename -- "$file")
	sampleID="${basefile%%.*}"
	echo $sampleID	$basefile	
	echo "$sampleID	$basefile" >> $pedfile	
done



