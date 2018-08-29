#!/bin/bash
#PBS -N GetGNOMAD_VNX
#PBS -V
#PBS -o /mnt/causes-vnx1/Databases/GetGNOMAD.o
#PBS -e /mnt/causes-vnx1/Databases/GetGNOMAD.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=128gb
## Set the max walltime for the job
#PBS -l walltime=200:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=32

source /opt/tools/hpcenv.sh


# Set your working directory
WORKING_DIR=/mnt/causes-vnx1/Databases/GNOMAD/RAW/
cd $WORKING_DIR
#
### 1) Get the data from online
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr1.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr2.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr3.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr4.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr5.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr6.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr7.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr8.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr9.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr10.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr11.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr12.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr13.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr14.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr15.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr16.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr17.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr18.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr19.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr20.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr21.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr22.vcf.bgz
#nohup wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chrX.vcf.bgz
#
## 2) normalize & decompose VCFs with vt
# Problems while using this ref, using a different one instead
#REF='/scratch/richmonp/GENOME/ucsc.hg19.fasta'
REF=/mnt/causes-data01/data/GENOMES/GSC/GRCh37-lite.fa
for i in {1..22}
do 
	VCFGZ='gnomad.genomes.r2.0.2.sites.chr'$i'.vcf.bgz'
	/opt/tools/tabix/tabix $VCFGZ
	NORMVCF='gnomad.genomes.r2.0.2.sites.chr'$i'.norm.vcf.bgz'
	
	zless $VCFGZ  \
	   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	   | /opt/tools/vt/vt decompose -s - \
	   | /opt/tools/vt/vt normalize -r $REF - \
   | /opt/tools/tabix/bgzip -c > $NORMVCF 
	/opt/tools/tabix/tabix -p vcf $NORMVCF
done

VCFGZ='gnomad.genomes.r2.0.2.sites.chrX.vcf.bgz'
NORMVCF='gnomad.genomes.r2.0.2.sites.chrX.norm.vcf.bgz'

zless $VCFGZ  \
           | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
           | /opt/tools/vt/vt decompose -s - \
           | /opt/tools/vt/vt normalize -r $REF - \
           | /opt/tools/tabix/bgzip -c > $NORMVCF
/opt/tools/tabix/tabix -p vcf $NORMVCF

# 3) Concatenate Normalized files
/opt/tools/bcftools/bin/bcftools concat gnomad.genomes.r2.0.2.sites.chr1.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr2.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr3.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr4.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr5.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr6.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr7.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr8.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr9.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr10.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr11.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr12.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr13.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr14.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr15.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr16.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr17.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr18.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr19.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr20.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr21.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr22.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chrX.norm.vcf.bgz | /opt/tools/tabix/bgzip -c > gnomad.genomes.r2.0.2.sites.wholeGenome.norm.vcf.bgz

/opt/tools/tabix/tabix gnomad.genomes.r2.0.2.sites.wholeGenome.norm.vcf.bgz

