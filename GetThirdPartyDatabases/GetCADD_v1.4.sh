#!/bin/bash
#PBS -N GetGNOMAD_VNX
#PBS -V
#PBS -o /mnt/causes-vnx1/DATABASES/GetGNOMAD_v2.1.o
#PBS -e /mnt/causes-vnx1/DATABASES/GetGNOMAD_v2.1.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=128gb
## Set the max walltime for the job
#PBS -l walltime=200:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=24

source /opt/tools/hpcenv.sh

# Set your working directory
WORKING_DIR=/mnt/causes-vnx1/DATABASES/CADD/v1.4/
cd $WORKING_DIR

### 1) Get the data from online
wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz.tbi
wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/gnomad.genomes.r2.0.1.sites.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/gnomad.genomes.r2.0.1.sites.tsv.gz.tbi

