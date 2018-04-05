#!/bin/bash
#SBATCH --account=rrg-wyeth
## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL
## CPU Usage
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --time=2-0:00
#SBATCH --nodes=1
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


# Set your working directory
WORKING_DIR=/home/richmonp/scratch/DATABASES/GNOMAD/RAW/
cd $WORKING_DIR

## 1) Get the data from online
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

# 2) normalize & decompose VCFs with vt
# Problems while using this ref, using a different one instead
#REF='/scratch/richmonp/GENOME/ucsc.hg19.fasta'
REF=/home/richmonp/project/GENOME/human_g1k_v37.fasta
for i in {1..22}
do 
	VCFGZ='gnomad.genomes.r2.0.2.sites.chr'$i'.vcf.bgz'
	tabix $VCFGZ
	NORMVCF='gnomad.genomes.r2.0.2.sites.chr'$i'.norm.vcf.bgz'
	
	zless $VCFGZ  \
	   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	   | vt decompose -s - \
	   | vt normalize -r $REF - \
	   | bgzip -c > $NORMVCF 
	tabix -p vcf $NORMVCF
done

VCFGZ='gnomad.genomes.r2.0.2.sites.chrX.vcf.bgz'
NORMVCF='gnomad.genomes.r2.0.2.sites.chrX.norm.vcf.bgz'

zless $VCFGZ  \
           | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
           | vt decompose -s - \
           | vt normalize -r $REF - \
           | bgzip -c > $NORMVCF
tabix -p vcf $NORMVCF

# 3) Concatenate Normalized files
bcftools concat gnomad.genomes.r2.0.2.sites.chr1.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr2.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr3.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr4.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr5.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr6.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr7.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr8.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr9.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr10.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr11.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr12.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr13.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr14.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr15.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr16.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr17.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr18.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr19.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr20.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr21.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chr22.norm.vcf.bgz gnomad.genomes.r2.0.2.sites.chrX.norm.vcf.bgz | bgzip -c > gnomad.genomes.r2.0.2.sites.wholeGenome.norm.vcf.bgz

tabix gnomad.genomes.r2.0.2.sites.wholeGenome.norm.vcf.bgz

