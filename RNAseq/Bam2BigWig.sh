#!/bin/bash

#SBATCH --partition=defq

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=160G
# Try to keep RAM=cores*8, since there is roughly 4G RAM/core

#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


# Bam2matrix for RNAseq
source /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda activate DeepTools


GFF3=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/Homo_sapiens.GRCh38.100.chr.gff3.gz
GTF=/mnt/common/DATABASES/REFERENCES/GRCh38/GENOME/GENCODE/gencode.v38.annotation.gtf
Working_Dir=/mnt/scratch/Precision/Hub/PROCESS/DH3132/
BAM1=F58387_4_lanes_dupsFlagged.bam
BAM2=F58388_4_lanes_dupsFlagged.bam
BAM3=F58389_4_lanes_dupsFlagged.bam
Family_ID=DH3132

cd $Working_Dir

bamCoverage -p 40 \
	-b $BAM1 \
	--normalizeUsing CPM \
	-o F58387.bw

bamCoverage -p 40 \
	--normalizeUsing CPM \
	-b $BAM2 \
	-o F58388.bw

bamCoverage -p 40 \
	--normalizeUsing CPM \
	-b $BAM3 \
	-o F58389.bw

