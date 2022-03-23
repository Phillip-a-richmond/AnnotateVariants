#!/bin/bash

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


##########
# Set up #
##########


# open up scratch
sudo chmod ugo=rwx -R /scratch/

# Running for this sample, Ideally would automate this from batch script to pull from file
Sample_ID=NA21143

## Set working space
Working_Dir=/scratch/$Sample_ID
mkdir -p $Working_Dir
cd $Working_Dir

# Where we pull the reads to, this is on the scratch space (onboard NVMe) 
CRAM_Dir=$Working_Dir/${Sample_ID}_CRAM
mkdir -p $CRAM_Dir

Sample_CRAM=${Sample_ID}.alt_bwamem_GRCh38DH.20150718.GIH.low_coverage.cram
Sample_VCF=${Sample_ID}.vcf.gz
Sample_GVCF=${Sample_ID}.gvcf.gz

Output_Dir=$Working_Dir/Variants/
mkdir -p $Output_Dir

### For final output
Final_Dir=/shared/DeepVariantOutput/

#### GRCh38 #### 
echo "GRCh38 genome"
Genome=GRCh38
Seq_Type=WGS
Fasta_Dir=/shared/Genomes/
Fasta_File=GRCh38_full_analysis_set_plus_decoy_hla.fa


# Pull data from s3
## Echo out a dry run
echo "Pulling data"
sudo aws s3 cp --no-sign-request --recursive --dryrun s3://1000genomes/data/${Sample_ID}/alignment/ $CRAM_Dir

# Run it
sudo aws s3 cp --no-sign-request --recursive s3://1000genomes/data/${Sample_ID}/alignment/ $CRAM_Dir

# Test you have it
ls $CRAM_Dir/$Sample_CRAM

# Get Docker
BIN_VERSION="1.1.0"
sudo apt -y update
sudo apt-get -y install docker.io
sudo docker pull google/deepvariant:"${BIN_VERSION}"

# Call variants
#sudo docker run \
#	-v "${CRAM_Dir}":"/cramdir" \
#	-v "${Fasta_Dir}":"/genomedir" \
#	-v "${Output_Dir}":"/output" \
#	google/deepvariant:"${BIN_VERSION}" \
#  /opt/deepvariant/bin/run_deepvariant \
#  --model_type=${Seq_Type} \
#  --ref="/genomedir/$Fasta_File" \
#	  --intermediate_results_dir="/output/intermediate_results_dir" \
#  --reads="/cramdir/$Sample_CRAM" \
#  --output_vcf="/output/${Sample_VCF}" \
#  --output_gvcf="/output/${Sample_GVCF}" \
#  --num_shards=$SLURM_CPUS_PER_TASK 

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


# Copy data back to final directory
cp $Output_Dir/${Sample_GVCF} $Final_Dir
cp $Output_Dir/${Sample_VCF} $Final_Dir

cp $Sample_ID.bed.gz $Final_Dir

exit

