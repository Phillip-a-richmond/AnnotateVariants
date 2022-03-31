#!/bin/bash


## CPU Usage
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error
# #SBATCH --array=2-21%20

##########
# Set up #
##########


# open up scratch
sudo chmod ugo=rwx -R /scratch/
sudo chmod ugo=rwx -R /shared/

##########
# Part 1 #
##########

# Here we will arrange which trio sample we are working on

# Running for this trio stored in a ped file. All the ped files are stored in that directory, which we will loop through in the next script
# Setting this for dev before running the array
SLURM_ARRAY_TASK_ID=1

# Found this online. real useful
# https://stackoverflow.com/questions/21668471/bash-script-create-array-of-all-files-in-a-directory
shopt -s nullglob

# Pulling from our list of ped files here
pedfiles=(/shared/AnnotateVariants/Cloud/1kG_Data/*ped)
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
Working_Dir=/shared/$childid
mkdir -p $Working_Dir
cd $Working_Dir

##########
# Part 2 #
##########

# Download the data for father, mother, child from s3, and rename it

# download father data
aws s3 cp --dryrun --recursive --no-sign-request s3://1000genomes/data/$fatherid/alignment/ .
aws s3 cp --recursive --no-sign-request s3://1000genomes/data/$fatherid/alignment/ .
# rename it
mv $fatherid*cram $fatherid.cram
mv $fatherid*cram.crai $fatherid.cram.crai

# download mother data
aws s3 cp --dryrun --recursive --no-sign-request s3://1000genomes/data/$motherid/alignment/ .
aws s3 cp --recursive --no-sign-request s3://1000genomes/data/$motherid/alignment/ .
# rename it
mv $motherid*cram $motherid.cram
mv $motherid*cram.crai $motherid.cram.crai

# download child data
aws s3 cp --recursive --no-sign-request --dryrun  s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ ./ --exclude "*" --include "*$childid*"
aws s3 cp --recursive --no-sign-request  s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ ./ --exclude "*" --include "*$childid*"
mv ./ERR*/$childid*cram $childid.cram
mv ./ERR*/$childid*cram.crai $childid.cram.crai


# Where we pull the reads to, this is on the scratch space (onboard NVMe) 
CRAM_Dir=$Working_Dir

Output_Dir=$Working_Dir/Variants/
mkdir -p $Output_Dir

### For final output
Final_Dir=/shared/SVOutput/
mkdir -p $Final_Dir

#### GRCh38 #### 
echo "GRCh38 genome"
Genome=GRCh38
Seq_Type=WGS
Fasta_Dir=/shared/AnnotateVariants/Cloud/Genomes/
Fasta_File=GRCh38_full_analysis_set_plus_decoy_hla.fa

##########
# Part 3 #
##########

# here I will Manta joint for the trio

## Activate it
source /shared/AnnotateVariants/Cloud/miniconda3/etc/profile.d/conda.sh
conda activate /shared/AnnotateVariants/Cloud/miniconda3/envs/Mamba/envs/SeqTools


## Step 1 - Make Config Script
configManta.py \
        --referenceFasta=$Fasta_Dir/$Fasta_File \
        --runDir=$Working_Dir \
        --bam $CRAM_Dir/$childid.cram \
        --bam $CRAM_Dir/$motherid.cram \
        --bam $CRAM_Dir/$fatherid.cram

## Step 2 - Execute config script
cd $Working_Dir
./runWorkflow.py \
        -j 32 \
        -g 64


##########
# Part 4 #
##########

# Get Docker
sudo apt -y update
sudo apt-get -y install docker.io
sudo docker pull brentp/smoove

# get the smoove.sh script
cp /shared/AnnotateVariants/Cloud/smoove.sh $Working_Dir

# Father Call variants
sudo docker run \
        -v "${CRAM_Dir}":"/cramdir" \
        -v "${Fasta_Dir}":"/genomedir" \
        -v "${Output_Dir}":"/output" \
        -it brentp/smoove \
        bash /bamdir/smoove.sh 32 $childid $Fasta_File \


##########
# Part 5 #
##########


# Copy Variant data back to /shared/
## Smoove
cp $Output_Dir/results-smoove/${childid}-smoove.genotyped.vcf.gz $Final_Dir

## Manta
cp $Output_Dir/results/variants/diploidSV.vcf.gz $Final_Dir/${childid}_Manta_diploidSV.vcf.gz


exit

##########
# Part 6 #
##########

# Get discordant reads in bed file

## Get tool
wget -O excord https://github.com/brentp/excord/releases/download/v0.2.2/excord_linux64
chmod +x ./excord

## Activate miniconda environment
source /shared/miniconda3/etc/profile.d/conda.sh 
conda activate /shared/miniconda3/envs/Mamba/envs/SeqTools

# Make a bam file
samtools view -@ 8 -b $CRAM_Dir/$Sample_CRAM -o $CRAM_Dir/$Sample_ID.bam

samtools view -b $CRAM_Dir/$Sample_ID.bam |
        ./excord \
        --discordantdistance 500 \
        --fasta $Fasta_Dir/$Fasta_File \
        /dev/stdin \
        | LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n -k3,3n \
        | bgzip -c $Output_Dir/$Sample_ID.bed.gz


cp $Sample_ID.bed.gz $Final_Dir

exit
