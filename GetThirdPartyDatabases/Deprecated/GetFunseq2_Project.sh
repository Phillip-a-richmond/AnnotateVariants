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
cd /project/projects/def-wyeth/DATABASES/

wget http://archive.gersteinlab.org/funseq2.1.2/hg19_NCscore_funseq216.tsv.bgz
md5sum hg19_NCscore_funseq216.tsv.bgz > hg19_NCscore_funseq216.tsv.bgz.md5sum
wget http://archive.gersteinlab.org/funseq2.1.2/hg19_NCscore_funseq216.tsv.bgz.tbi
md5sum hg19_NCscore_funseq216.tsv.bgz.tbi > hg19_NCscore_funseq216.tsv.bgz.tbi.md5sum

