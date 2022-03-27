#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error
#SBATCH --array=2-21%20

##########
# Set up #
##########


# open up scratch
sudo chmod ugo=rwx -R /scratch/


##########
# Part 1 #
##########

# Here we will arrange which trio sample we are working on

# Running for this trio stored in a ped file. All the ped files are stored in that directory, which we will loop through in the next script
# Setting this for dev before running the array
#SLURM_ARRAY_TASK_ID=1

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
Working_Dir=/scratch/$childid
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
Final_Dir=/shared/DeepVariantOutput/

#### GRCh38 #### 
echo "GRCh38 genome"
Genome=GRCh38
Seq_Type=WGS
Fasta_Dir=/shared/AnnotateVariants/Cloud/Genomes/
Fasta_File=GRCh38_full_analysis_set_plus_decoy_hla.fa

##########
# Part 3 #
##########

# here I will get docker, then run deepvariant separately for each sample

# Get Docker
BIN_VERSION="1.3.0"
sudo apt -y update
sudo apt-get -y install docker.io
sudo docker pull google/deepvariant:"${BIN_VERSION}"

# Father Call variants
sudo docker run \
	-v "${CRAM_Dir}":"/cramdir" \
	-v "${Fasta_Dir}":"/genomedir" \
	-v "${Output_Dir}":"/output" \
	google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=${Seq_Type} \
  --ref="/genomedir/$Fasta_File" \
  --intermediate_results_dir="/output/intermediate_results_dir" \
  --reads="/cramdir/$fatherid.cram" \
  --output_vcf="/output/$fatherid.deepvariant.$BIN_VERSION.vcf.gz" \
  --output_gvcf="/output/$fatherid.deepvariant.$BIN_VERSION.gvcf.gz" \
  --num_shards=16 

# Mother Call variants
sudo docker run \
        -v "${CRAM_Dir}":"/cramdir" \
        -v "${Fasta_Dir}":"/genomedir" \
        -v "${Output_Dir}":"/output" \
        google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=${Seq_Type} \
  --ref="/genomedir/$Fasta_File" \
  --intermediate_results_dir="/output/intermediate_results_dir" \
  --reads="/cramdir/$motherid.cram" \
  --output_vcf="/output/$motherid.deepvariant.$BIN_VERSION.vcf.gz" \
  --output_gvcf="/output/$motherid.deepvariant.$BIN_VERSION.gvcf.gz" \
  --num_shards=16    

# Child Call variants
sudo docker run \
        -v "${CRAM_Dir}":"/cramdir" \
        -v "${Fasta_Dir}":"/genomedir" \
        -v "${Output_Dir}":"/output" \
        google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=${Seq_Type} \
  --ref="/genomedir/$Fasta_File" \
  --intermediate_results_dir="/output/intermediate_results_dir" \
  --reads="/cramdir/$childid.cram" \
  --output_vcf="/output/$childid.deepvariant.$BIN_VERSION.vcf.gz" \
  --output_gvcf="/output/$childid.deepvariant.$BIN_VERSION.gvcf.gz" \
  --num_shards=16


##########
# Part 4 #
##########

# Run deepTrio
sudo docker pull google/deepvariant:deeptrio-"${BIN_VERSION}"

sudo docker run \
  -v "${CRAM_Dir}":"/cramdir" \
  -v "${Fasta_Dir}":"/genomedir" \
  -v "${Output_Dir}":"/output" \
  google/deepvariant:deeptrio-"${BIN_VERSION}" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type=WGS \
  --intermediate_results_dir="/output/intermediate_results_dir" \
  --ref="/genomedir/$Fasta_File" \
  --sample_name_child "$childid" \
  --sample_name_parent1 "$fatherid" \
  --sample_name_parent2 "$motherid" \
  --reads_child=/cramdir/$childid.cram \
  --reads_parent1=/cramdir/$fatherid.cram \
  --reads_parent2=/cramdir/$motherid.cram \
  --output_vcf_child /output/$childid.deeptrio.$BIN_VERSION.vcf.gz \
  --output_vcf_parent1 /output/$fatherid.deeptrio.$BIN_VERSION.vcf.gz \
  --output_vcf_parent2 /output/$motherid.deeptrio.$BIN_VERSION.vcf.gz \
  --output_gvcf_child /output/$childid.deeptrio.$BIN_VERSION.gvcf.gz \
  --output_gvcf_parent1 /output/$fatherid.deeptrio.$BIN_VERSION.gvcf.gz \
  --output_gvcf_parent2 /output/$motherid.deeptrio.$BIN_VERSION.gvcf.gz \
  --num_shards 16  


##########
# Part 5 #
##########

# Copy Variant data back to /shared/
cp $Output_Dir/*vcf.gz $Final_Dir
cp $Output_Dir/*gvcf.gz $Final_Dir
cp $Output_Dir/*tbi $Final_Dir
cp $Output_Dir/*html $Final_Dir


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

