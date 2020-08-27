# This script will fetch the fasta files 
conda activate STR_environment

# Where do you want them
GENOME_DIR=/project/st-wasserww-1/GENOME/
mkdir -p $GENOME_DIR
cd $GENOME_DIR


###############################################################################
# Step 2 - Get Reference genomes and index them. #
##################################################

# Get Reference Genome Files (these match the files in str_bams) 
## GRCh37
mkdir -p GRCh37
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa.fai

bwa index GRCh37-lite.fa

## GRCh38
mkdir -p GRCh38
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz > GRCh38-lite.fa

samtools faidx GRCh38-lite.fa
bwa index GRCh38-lite.fa


