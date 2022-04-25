#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=160G
#SBATCH --cpus-per-task=20
#SBATCH --time=8:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

# Load singularity
module load singularity

Stix_SIF=/mnt/common/Precision/STIX/stix.sif

WorkingDir=/mnt/common/OPEN_DATA/1kG_Trio/
DataDir=$WorkingDir/Excord_Output/
SVDir=$WorkingDir/MantaSmoove/
IndexDir=$WorkingDir/STIX_Index_Try2/
PedFile=$DataDir/AllSamples.ped

# Initialize stix database with giggle index and bed files
singularity exec \
	-B $WorkingDir \
	-B $IndexDir \
	-B $DataDir \
	$Stix_SIF \
	stix \
	-i $IndexDir \
	-p $PedFile \
	-d $DataDir/AllSamples.stix.db \
	-c 2


