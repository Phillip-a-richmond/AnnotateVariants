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
WORKING_DIR=/mnt/causes-vnx1/DATABASES/GNOMAD/v2.1/
cd $WORKING_DIR
#
### 1) Get the data from online
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz.tbi

## 2) normalize & decompose VCFs with vt
# Problems while using this ref, using a different one instead
REF=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa
VCFGZ=gnomad.genomes.r2.1.sites.vcf.bgz
NORMVCF=gnomad.genomes.r2.1.sites.norm.vcf.bgz
zless $VCFGZ  \
   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
   | /opt/tools/vt/vt decompose -s - \
   | /opt/tools/vt/vt normalize -r $REF - \
   | /opt/tools/tabix/bgzip -c > $NORMVCF 
/opt/tools/tabix/tabix -p vcf $NORMVCF

# Repeat for exome
#1 get data 
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz.tbi

# 2 normalize and decompose
REF=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa
VCFGZ=gnomad.exomes.r2.1.sites.vcf.bgz
NORMVCF=gnomad.exomes.r2.1.sites.norm.vcf.bgz
zless $VCFGZ  \
   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
   | /opt/tools/vt/vt decompose -s - \
   | /opt/tools/vt/vt normalize -r $REF - \
   | /opt/tools/tabix/bgzip -c > $NORMVCF 
/opt/tools/tabix/tabix -p vcf $NORMVCF


