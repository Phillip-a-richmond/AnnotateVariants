#!/bin/bash
#SBATCH --account=rrg-wyeth
## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL
## CPU Usage
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --time=4-0:00
#SBATCH --nodes=1
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error
# Load Modules
cd /home/richmonp/scratch/DATABASES/

nohup wget http://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_noncoding.vcf.gz
md5sum fathmm_xf_noncoding.vcf.gz > fathmm_xf_noncoding.vcf.gz.md5sum
nohup wget http://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_noncoding.vcf.gz.tbi

nohup wget http://fathmm.biocompute.org.uk/database/fathmm-MKL_Current_zerobased.tab.gz
md5sum fathmm-MKL_Current_zerobased.tab.gz > fathmm-MKL_Current_zerobased.tab.gz.md5sum
tabix -p bed fathmm-MKL_Current_zerobased.tab.gz
