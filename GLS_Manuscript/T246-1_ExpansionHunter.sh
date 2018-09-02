#!/bin/bash
#PBS -N STR_Genotyping_GLS
#PBS -V
#PBS -o /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T246/T246-1.o
#PBS -e /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T246/T246-1.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=30gb
## Set the max walltime for the job
#PBS -l walltime=200:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=8
NSLOTS=$PBS_NUM_PPN
umask 0002
source /opt/tools/hpcenv.sh

SAMPLE_ID='TIDEX913_BWAmem'
WORKING_DIR='/mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T246/'
GENOME_FASTA='/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa'

#EH
/opt/tools/ExpansionHunter-v2.5.3-linux_x86_64/bin/ExpansionHunter --bam $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' \
	--repeat-specs /opt/tools/ExpansionHunter-v2.5.3-linux_x86_64/data/repeat-specs/grch37 \
	--ref-fasta $GENOME_FASTA \
	--log $WORKING_DIR${SAMPLE_ID}_ExpansionHunter.log \
	--json $WORKING_DIR${SAMPLE_ID}_ExpansionHunter.json \
	--vcf $WORKING_DIR${SAMPLE_ID}_ExpansionHunter.vcf \
	--sex male

#EH DN
/mnt/causes-vnx1/PIPELINES/ExpansionHunterDenovo-v0.6.2 --bam $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \
        --reference $GENOME_FASTA \
        --output $WORKING_DIR${SAMPLE_ID}_ExpansionHunterDenovo.json

# GangSTR GenomeWide
GangSTR \
        --bam $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \
        --ref $GENOME_FASTA \
        --genomewide \
        --out ${SAMPLE_ID}_GangSTR_Genomewide \
        --regions /mnt/causes-vnx1/DATABASES/GRCh37_ver8.bed	

# GangSTR local
GangSTR \
	--bam $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \
	--ref $GENOME_FASTA \
	--regions /mnt/causes-vnx1/DATABASES/GLS.bed --out ${SAMPLE_ID}_GLS_GangSTR 

GangSTR \
	--bam $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \
	--ref $GENOME_FASTA \
	--readlength 150 \
	--regions /mnt/causes-vnx1/DATABASES/GLS.bed --out ${SAMPLE_ID}_GLS_GangSTR_readlen150 \



