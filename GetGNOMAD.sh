#!/bin/bash
#SBATCH --account=rrg-wyeth
## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL
## CPU Usage
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --time=2-0:00
#SBATCH --nodes=1
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error



# Load Modules

cd /home/richmonp/scratch/DATABASES/




# Get the data from online
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.1.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.2.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.3.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.4.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.5.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.6.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.7.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.8.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.9.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.10.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.11.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.12.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.13.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.14.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.15.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.16.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.17.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.18.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.19.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.20.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.21.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.22.vcf.gz
nohup wget https://storage.googleapis.com/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.X.vcf.gz

# 1) normalize & decompose VCFs with vt

REF='/mnt/data/GENOMES/GSC/GRCh37-lite.fa'
for i in {1..22}
do 
	VCFGZ='gnomad.genomes.r2.0.1.sites.'$i'.vcf.gz'
	NORMVCF='gnomad.genomes.r2.0.1.sites.'$i'.norm.vcf.gz'
	
	zless $VCFGZ  \
	   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	   | vt decompose -s - \
	   | vt normalize -r $REF - \
	   | bgzip -c > $NORMVCF 
	tabix -p vcf $NORMVCF
done

VCFGZ='gnomad.genomes.r2.0.1.sites.X.vcf.gz'
NORMVCF='gnomad.genomes.r2.0.1.sites.X.norm.vcf.gz'

zless $VCFGZ  \
           | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
           | vt decompose -s - \
           | vt normalize -r $REF - \
           | bgzip -c > $NORMVCF
tabix -p vcf $NORMVCF

# 2) Concatenate Normalized files
vcf-concat gnomad.genomes.r2.0.1.sites.1.norm.vcf.gz gnomad.genomes.r2.0.1.sites.2.norm.vcf.gz gnomad.genomes.r2.0.1.sites.3.norm.vcf.gz gnomad.genomes.r2.0.1.sites.4.norm.vcf.gz gnomad.genomes.r2.0.1.sites.5.norm.vcf.gz gnomad.genomes.r2.0.1.sites.6.norm.vcf.gz gnomad.genomes.r2.0.1.sites.7.norm.vcf.gz gnomad.genomes.r2.0.1.sites.8.norm.vcf.gz gnomad.genomes.r2.0.1.sites.9.norm.vcf.gz gnomad.genomes.r2.0.1.sites.10.norm.vcf.gz gnomad.genomes.r2.0.1.sites.11.norm.vcf.gz gnomad.genomes.r2.0.1.sites.12.norm.vcf.gz gnomad.genomes.r2.0.1.sites.13.norm.vcf.gz gnomad.genomes.r2.0.1.sites.14.norm.vcf.gz gnomad.genomes.r2.0.1.sites.15.norm.vcf.gz gnomad.genomes.r2.0.1.sites.16.norm.vcf.gz gnomad.genomes.r2.0.1.sites.17.norm.vcf.gz gnomad.genomes.r2.0.1.sites.18.norm.vcf.gz gnomad.genomes.r2.0.1.sites.19.norm.vcf.gz gnomad.genomes.r2.0.1.sites.20.norm.vcf.gz gnomad.genomes.r2.0.1.sites.21.norm.vcf.gz gnomad.genomes.r2.0.1.sites.22.norm.vcf.gz gnomad.genomes.r2.0.1.sites.X.norm.vcf.gz | bgzip -c > gnomad.genomes.r2.0.1.sites.wholeGenome.norm.vcf.gz

