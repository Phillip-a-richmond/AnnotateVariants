# This script cleans up the emtpy excord files.


## Define directory for working
WorkingDir=/mnt/common/OPEN_DATA/1kG_Trio/
cd $WorkingDir

## Define directory for output
EmptyDir=$WorkingDir/FailedRuns/
mkdir -p $EmptyDir

# This one liner does what I want
# https://superuser.com/questions/644272/how-do-i-delete-all-files-smaller-than-a-certain-size-in-all-subfolders

# this one has find and move
# https://unix.stackexchange.com/questions/484791/bash-find-move-files-bigger-than-certain-size 

# Find the files with size less than 50k
find . -name "*.bed.gz" -type 'f' -size -50k 

# Find them and move them to this empty directory
find . -name "*.bed.gz" -type 'f' -size -50k -exec mv "{}" $EmptyDir \;


# For the remaining files, lets loop through them and get their filesizes

# Used for dev
# File=NA19203.bed.gz

# Found this online. real useful
# https://stackoverflow.com/questions/21668471/bash-script-create-array-of-all-files-in-a-directory
shopt -s nullglob
Files=($WorkingDir/*.bed.gz)

# Store record of the emtpies
EmptyLog=$WorkingDir/EmptyExcordRuns.txt
FullLog=$WorkingDir/FullExcordRuns.txt

for File in ${Files[@]}; 
do
	fileSize=`du -k "$File" | cut -f1`
	echo $File
	echo $fileSize

	## From https://unix.stackexchange.com/questions/6758/how-can-i-check-if-a-gzipped-file-is-empty
	#if LC_ALL=C gzip -l "$File" | awk 'NR==2 {exit($2!=0)}'; then
	#  echo "$File is empty"
	#  echo "$File" >> $EmptyLog
	#  mv $File $EmptyDir
	#else
	#  echo "$File is not empty"
	#  fileSize=$(stat -c%s "$File")
	#  echo "$File	$fileSize" >> $FullLog 
	#fi
done

