#!/bin/sh
## Author: Phillip Richmond (prichmond@bcchr.ca) 
## LICENSE: CC-BY-SA 4.0. https://creativecommons.org/licenses/by-sa/4.0/ 
#
#######################################
## Step 0 - get miniconda3 and docker #
#######################################
#
## This Gist assumes you have Miniconda. If no Miniconda, fix your life, do step 0.
  
function buildMiniConda3() 
{
        DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
        MINI_CONDA_INSTALL_DIR=$DIR/miniconda3

        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b -p $MINI_CONDA_INSTALL_DIR

        source ${MINI_CONDA_INSTALL_DIR}/etc/profile.d/conda.sh
        cd $DIR
        conda --help
        rm Miniconda3-latest-Linux-x86_64.sh
}

# Execute it, or don't if you can point to a conda
buildMiniConda3

##########################################
## Step 1 set your miniconda3, get tools #
##########################################
#
### If you ran setup above, use this
MINICONDA3=$PWD/miniconda3
#
### Otherwise, where is your miniconda3
##MINICONDA3=/mnt/common/Precision/Miniconda3/opt/miniconda3/
source $MINICONDA3/etc/profile.d/conda.sh
#
### Mamba is better, will help solve envs far faster
#conda create -y -c conda-forge -n Mamba \
#	mamba
#	
#conda activate Mamba
#
### Get your tools. First time I ran this it took awhile to solve, so kept fgbio separate for now.
#mamba create -y -c bioconda -n FGbio \
#	fgbio=1.5.1 
#
#### here I'm adding bwa,samtools, htslib, picard
#mamba create -y -c bioconda -n SeqTools \
#	samtools htslib bwa picard 
#
### Activate it to get our tools
### Note, put this inside a mamba env, so need to specify the env path
#conda activate $MINICONDA3/envs/Mamba/envs/FGbio
#fgbio --help
#
### Test this environment too
#conda activate $MINICONDA3/envs/Mamba/envs/SeqTools
#bwa --help
#samtools --help
#
#
########################
# Step 2: Get new data #
########################

## Get new Fasta
function GetFastas()
{
	mkdir $PWD/Genomes/
	cd $PWD/Genomes/
	wget -c -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz 
	gunzip GCA_000001405.29_GRCh38.p14_genomic.fna.gz 

	## Get assembly report
	wget -c -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_assembly_report.txt

	## Get new gff
	### Note: don't think these will be all that different, and mostly use VEP to manage this. Coordinate system doesn't change, so no real worries here.
	### Only do this if you have to
	#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.gff.gz
	
	# Get 1000Genomes too
	wget -c -q http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
}

GetFastas 

################################
# Step 3: Create a fasta index #
################################
## Activate tools
conda activate $MINICONDA3/envs/Mamba/envs/SeqTools

## Samtools faidx
cd $PWD/Genomes/
samtools faidx GCA_000001405.29_GRCh38.p14_genomic.fna

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa


###############################################
# Step 4: Create a seq dict we currently have #
###############################################

## Picard create sequence dictionary
picard CreateSequenceDictionary \
	R=GCA_000001405.29_GRCh38.p14_genomic.fna \
	O=GCA_000001405.29_GRCh38.p14_genomic.dict

####################################################
# Step 5: create the seq dict we want, and sort it #
####################################################
conda activate $MINICONDA3/envs/Mamba/envs/FGbio

## Collect alternate contig names from the assembly report
fgbio CollectAlternateContigNames \
	-i GCA_000001405.29_GRCh38.p14_assembly_report.txt \
	-o GCA_000001405.29_GRCh38.p14_genomic_alternate.dict \
	-p Sequence-Name \
	-a GenBank-Accn 
	
## Sort the sequence dictionary we just got from the report, to match our current (downloaded GRCh38.14) fasta dict order
fgbio SortSequenceDictionary \
	-i GCA_000001405.29_GRCh38.p14_genomic_alternate.dict \
	-d GCA_000001405.29_GRCh38.p14_genomic.dict \
	-o GCA_000001405.29_GRCh38.p14_genomic_sorted_alternate.dict

############################################
# Step 6: change the seqnames of the fasta #
############################################
fgbio UpdateFastaContigNames \
	-i GCA_000001405.29_GRCh38.p14_genomic.fna \
	-d GCA_000001405.29_GRCh38.p14_genomic_sorted_alternate.dict \
	-o GCA_000001405.29_GRCh38.p14_genomic_chromFixed.fa

