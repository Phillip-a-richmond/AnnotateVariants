#!/bin/bash

# Found this online. real useful
# https://stackoverflow.com/questions/21668471/bash-script-create-array-of-all-files-in-a-directory
shopt -s nullglob

# Set directory with outfiles
outFileDir=$PWD
cd $outFileDir

outFiles=($outFileDir/*.out)

echo $outFiles

for file in ${outFiles[@]}
do
	echo $file
	grep ".ped" $file
	grep "Runtime"  $file 

done


