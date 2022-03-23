#!/bin/bash

SLURM_ARRAY_TASK_ID=0
shopt -s nullglob
pedfiles=(/shared/AnnotateVariants/Cloud/1kG_Data/*ped)
pedfile=${pedfiles[$SLURM_ARRAY_TASK_ID]}
echo $pedfile
fatherid=$( awk 'NR==2 {print $1}' $pedfile)
motherid=$( awk 'NR==3 {print $1}' $pedfile)
childid=$( awk 'NR==4 {print $1}' $pedfile)

echo $fatherid
echo $motherid
echo $childid

# download father data
aws s3 cp --dryrun --recursive --no-sign-request s3://1000genomes/data/$fatherid/alignment/ .
#mv $fatherid*cram $fatherid.cram
#mv $fatherid*cram.crai $fatherid.cram.crai

# download mother data
aws s3 cp --dryrun --recursive --no-sign-request s3://1000genomes/data/$motherid/alignment/ .

# download child data
aws s3 cp --recursive --no-sign-request --dryrun  s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ ./ --exclude "*" --include "*$childid*"

exit

