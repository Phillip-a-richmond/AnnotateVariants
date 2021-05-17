#!/bin/bash
#SBATCH --account=rrg-wyeth
## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL
## CPU Usage
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH --time=2-0:00
#SBATCH --nodes=1
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

cd /home/richmonp/project/DATABASES/

wget http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg19/JASPAR2018_hg19_all_chr.bed.gz
md5sum JASPAR2018_hg19_all_chr.bed.gz > JASPAR2018_hg19_all_chr.bed.gz.md5sum

gunzip -c JASPAR2018_hg19_all_chr.bed.gz | bgzip  -c   >  JASPAR2018_hg19_all_chr.bed.bgz
tabix -p bed JASPAR2018_hg19_all_chr.bed.bgz

