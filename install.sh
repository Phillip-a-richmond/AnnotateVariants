#!/bin/bash

########################################################################
###
### Copyright (c) 2020, Wasserman lab
###
### FILE        install.sh
###
### DESCRIPTION This is the install.sh script for the project.
###             We use conda in order to install the application.
###
### Initial version @ Godfrain Jacques KOUNKOU
### Modified @ Phillip Richmond
########################################################################

# set -e

########
# Usage:
# bash install.sh -d </path/to/database/install/directory/>
# DESC This code block gets the database directory from the input options and stores it as DB_DIR
#########
usage()
{
    echo "usage: bash install.sh -d </path/to/database/install/directory/>"
    echo "example: $ bash install.sh -d /scratch/ex-ofornes-1/RICHMOND/ANNOTATE_VARIANTS/DATABASES/"
}

while getopts d: flag
do
    case "${flag}" in
        d) DB_DIR=${OPTARG};;
    esac
done

if [ -z "$DB_DIR" ]
then
    usage
    exit
fi

echo "Databases will be installed in: $DB_DIR";

##########################################
# DESC This function will find all tools needed for AnnotateVariants
# ARGS This function doesnt require any arguments
# RSLT Downloads all necessary tools (except java jdk and other tools listed below) into opt/AnnotateVariants
##########################################
function buildAnnotateVariants() 
{
	DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
	OPT_DIR=${DIR}/opt
	mkdir -p ${OPT_DIR}
	pushd ${OPT_DIR}
	MINI_CONDA_INSTALL_DIR=$OPT_DIR/miniconda3

	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p $MINI_CONDA_INSTALL_DIR 

	source ${MINI_CONDA_INSTALL_DIR}/etc/profile.d/conda.sh
	conda env create --prefix ${OPT_DIR}/AnnotateVariantsEnvironment  -f $DIR/AnnotateVariants_CondaEnv.yml 
	cd $DIR
	conda activate opt/AnnotateVariantsEnvironment
}


function createDatabaseDirectory()
{
	DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
	OPT_DIR=${DIR}/opt
	cd $DIR
	
	# We need bwa and samtools, so load the environment we created above
	MINI_CONDA_INSTALL_DIR=$OPT_DIR/miniconda3
	source ${MINI_CONDA_INSTALL_DIR}/etc/profile.d/conda.sh
	conda activate opt/AnnotateVariantsEnvironment

	mkdir -p $DB_DIR

	# Get Reference Genome Files 
	## GRCh37
	mkdir -p $DB_DIR/GRCh37/
	cd $DB_DIR/GRCh37/
	wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa
	wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa.fai

	bwa index GRCh37-lite.fa

	## GRCh38
	mkdir -p $DB_DIR/GRCh38/
	wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
	gunzip -c Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz > GRCh38-lite.fa

	samtools faidx GRCh38-lite.fa
	bwa index GRCh38-lite.fa

}


function getDatabases()
{
	echo "not yet"
}

buildAnnotateVariants
#createDatabaseDirectory
