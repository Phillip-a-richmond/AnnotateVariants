#!/bin/bash

## CPU Usage
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error
#SBATCH --array=101-500%100

##########
# Set up #
##########

# Where is AnnotateVariants
AnnotateVariantsDir=/shared/AnnotateVariants/

# For final output
Final_Dir=/shared/SVOutput/
sudo chmod ugo=rwx -R /shared/SVOutput/
mkdir -p $Final_Dir

# For active processing
ScratchDir=/scratch/


# Set a log file
# Timing from here to end of script
FullRunStart=`date +%s`
currentTime=`date`
echo "Timestamp. Starting analysis: $currentTime"
echo "In seconds: $FullRunStart"

##########
# Part 1 #
##########

# Here we will arrange which trio sample we are working on

# Running for this trio stored in a ped file. All the ped files are stored in that directory, which we will loop through in the next script
# Setting this for dev before running the array
# When running array, make sure you comment this out.
# SLURM_ARRAY_TASK_ID=0

# Found this online. real useful
# https://stackoverflow.com/questions/21668471/bash-script-create-array-of-all-files-in-a-directory
shopt -s nullglob

# Pulling from our list of ped files here
pedfiles=($AnnotateVariantsDir/Cloud/1kG_Data/*ped)
pedfile=${pedfiles[$SLURM_ARRAY_TASK_ID]}
echo $pedfile

# Use awk to get the sample ids. I know that the father is line 2 (NR==2) and is the first column ($1)
# mother is line 3, child is line 4
fatherid=$( awk 'NR==2 {print $1}' $pedfile)
motherid=$( awk 'NR==3 {print $1}' $pedfile)
childid=$( awk 'NR==4 {print $1}' $pedfile)

# Check the ids
echo $fatherid
echo $motherid
echo $childid

### Set working space
# move this back to scratch in production
Working_Dir=$ScratchDir/$childid

# open up scratch, only needed if running w/ onboard NVMe
sudo chmod ugo=rwx -R $ScratchDir

mkdir -p $Working_Dir
cd $Working_Dir


##########
# Part 2 #
##########

# Download the data for father, mother, child from s3, and rename it
Start=`date +%s`
currentTime=`date`
echo "Timestamp. Step2-Start: Downloading data, $currentTime"
echo "In seconds: $Start"


# download father data
aws s3 cp --dryrun --recursive --no-sign-request s3://1000genomes/data/$fatherid/alignment/ .
aws s3 cp --quiet --recursive --no-sign-request s3://1000genomes/data/$fatherid/alignment/ .
# rename it
mv $fatherid*cram $fatherid.cram
mv $fatherid*cram.crai $fatherid.cram.crai

# download mother data
aws s3 cp --dryrun --recursive --no-sign-request s3://1000genomes/data/$motherid/alignment/ .
aws s3 cp --quiet --recursive --no-sign-request s3://1000genomes/data/$motherid/alignment/ .
# rename it
mv $motherid*cram $motherid.cram
mv $motherid*cram.crai $motherid.cram.crai

# download child data
aws s3 cp --recursive --no-sign-request --dryrun  s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ ./ --exclude "*" --include "*$childid*"
aws s3 cp --quiet --recursive --no-sign-request  s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ ./ --exclude "*" --include "*$childid*"
mv ./ERR*/$childid*cram $childid.cram
mv ./ERR*/$childid*cram.crai $childid.cram.crai


# Where we pull the reads to, this is on the scratch space (onboard NVMe) 
CRAM_Dir=$Working_Dir

#### GRCh38 #### 
echo "GRCh38 genome"
Genome=GRCh38
Seq_Type=WGS
Fasta_Dir=$AnnotateVariantsDir/Cloud/Genomes/
Fasta_File=GRCh38_full_analysis_set_plus_decoy_hla.fa

End=`date +%s`
runtime=$((End-Start))
currentTime=`date`
echo "Timestamp. Step2-End: Downloading data, $currentTime"
echo "In seconds: $End"
echo "Step2 Runtime: $runtime"


##########
# Part 3 #
##########


Start=`date +%s`
currentTime=`date`
echo "Timestamp. Step3-Start child: Excord,  $currentTime"
echo "In seconds: $Start"

# Get discordant reads in bed file

## Activate miniconda environment
source $AnnotateVariantsDir/Cloud/miniconda3/etc/profile.d/conda.sh
conda activate $AnnotateVariantsDir/Cloud/miniconda3/envs/Mamba/envs/SeqTools
EXCORD=$AnnotateVariantsDir/Cloud/excord

## Child
Sample_ID=$childid
Sample_CRAM=$Sample_ID.cram

samtools view -b -u -T $Fasta_Dir/$Fasta_File $CRAM_Dir/$Sample_ID.cram | \
        $EXCORD \
        --discordantdistance 500 \
        --fasta $Fasta_Dir/$Fasta_File \
        /dev/stdin \
        | LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n -k3,3n \
        | bgzip -c > $CRAM_Dir/$Sample_ID.bed.gz

cp $CRAM_Dir/$Sample_ID.bed.gz $Final_Dir


End=`date +%s`
runtime=$((End-Start))
currentTime=`date`
echo "Timestamp. Step3-End child: Running excord,  $currentTime"
echo "In seconds: $End"
echo "Step3 Runtime: $runtime"


## Mother
Start=`date +%s`
currentTime=`date`
echo "Timestamp. Step3-Start mother: Excord,  $currentTime"
echo "In seconds: $Start"

Sample_ID=$motherid
Sample_CRAM=$Sample_ID.cram

samtools view -b -u -T $Fasta_Dir/$Fasta_File $CRAM_Dir/$Sample_ID.cram | \
        $EXCORD \
        --discordantdistance 500 \
        --fasta $Fasta_Dir/$Fasta_File \
        /dev/stdin \
        | LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n -k3,3n \
        | bgzip -c > $CRAM_Dir/$Sample_ID.bed.gz

cp $CRAM_Dir/$Sample_ID.bed.gz $Final_Dir

End=`date +%s`
runtime=$((End-Start))
currentTime=`date`
echo "Timestamp. Step3-End mother: Running excord,  $currentTime"
echo "In seconds: $End"
echo "Step3 Runtime: $runtime"

## Father
Start=`date +%s`
currentTime=`date`
echo "Timestamp. Step3-Start father: Excord,  $currentTime"
echo "In seconds: $Start"

Sample_ID=$fatherid
Sample_CRAM=$Sample_ID.cram

samtools view -b -u -T $Fasta_Dir/$Fasta_File $CRAM_Dir/$Sample_ID.cram | \
        $EXCORD \
        --discordantdistance 500 \
        --fasta $Fasta_Dir/$Fasta_File \
        /dev/stdin \
        | LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n -k3,3n \
        | bgzip -c > $CRAM_Dir/$Sample_ID.bed.gz

cp $CRAM_Dir/$Sample_ID.bed.gz $Final_Dir


End=`date +%s`
runtime=$((End-Start))
currentTime=`date`
echo "Timestamp. Step3-End father: Running excord,  $currentTime"
echo "In seconds: $End"
echo "Step3 Runtime: $runtime"

FullRunEnd=$End
currentTime=`date`
FullRuntime=$((FullRunEnd-FullRunStart))
echo "Timestamp. Final. $currentTime"
echo "In seconds: $FullRunEnd"
echo "Full Runtime: $FullRuntime"


exit

