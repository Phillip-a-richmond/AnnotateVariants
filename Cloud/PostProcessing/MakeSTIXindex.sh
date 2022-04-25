#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=160G
#SBATCH --cpus-per-task=20
#SBATCH --time=240:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

# Load singularity
module load singularity

Stix_SIF=/mnt/common/Precision/STIX/stix.sif

WorkingDir=/mnt/common/OPEN_DATA/1kG_Trio/
DataDir=$WorkingDir/Excord_Output/
IndexDir=${WorkingDir}STIX_Index
rm -rf $IndexDir
mkdir -p $IndexDir

# Make a giggle index bash script. This is necessary because singularity won't evaluate the blob format of $IndexDir/*gz
IndexScript=$WorkingDir/IndexGiggle.sh
# clean up if 1 exists already
rm $IndexScript

# here I'm making this giggle script that will be passed to singularity below
echo cd $DataDir > $IndexScript
echo giggle index -s -f -o $IndexDir -i \"$DataDir*.bed.gz\" >> $IndexScript

# Run giggle index
singularity exec \
	-B $WorkingDir \
	-B $IndexDir \
	-B $DataDir \
	$Stix_SIF \
	sh $IndexScript 

