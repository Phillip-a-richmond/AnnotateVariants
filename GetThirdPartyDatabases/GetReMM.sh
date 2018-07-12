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

wget http://remm.visze.de/files/ReMM.v0.3.1.tsv.gz
md5sum ReMM.v0.3.1.tsv.gz > ReMM.v0.3.1.tsv.gz.md5sum
tabix -s 1 -b 2 ReMM.v0.3.1.tsv.gz
