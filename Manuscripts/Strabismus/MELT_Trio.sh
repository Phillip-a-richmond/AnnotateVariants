# Check to load environment
source /opt/tools/hpcenv.sh

WORKING_DIR=/mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T308/
ANALYSIS_DIR=${WORKING_DIR}MEI/
BAM1=12D3695_dedup.bam
BAM2=12D3698_dedup.bam
BAM3=17D8290_dedup.bam
MELT_DIR=/opt/tools/MELTv2.1.5/
MEI_LIST=${MELT_DIR}/me_refs/1KGP_Hg19/mei_list.txt
GENOME_FASTA=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa
GENE_ANNO=/opt/tools/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed

# Step 1 - Preprocess
java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
        -bamfile $WORKING_DIR$BAM1 \
	-h $GENOME_FASTA

java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
        -bamfile $WORKING_DIR$BAM2 \
	-h $GENOME_FASTA

java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
        -bamfile $WORKING_DIR$BAM3 \
	-h $GENOME_FASTA


# Step 2 - Individual Analysis
java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
	-w $ANALYSIS_DIR \
	-t $MEI_LIST \
	-c 30 \
	-bamfile $WORKING_DIR$BAM1 \
	-h $GENOME_FASTA

java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
        -w $ANALYSIS_DIR \
        -t $MEI_LIST \
        -c 30 \
        -bamfile $WORKING_DIR$BAM2 \
        -h $GENOME_FASTA

java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
        -w $ANALYSIS_DIR \
        -t $MEI_LIST \
        -c 30 \
        -bamfile $WORKING_DIR$BAM3 \
        -h $GENOME_FASTA


# Step 3 Group Discovery
java -Xmx4G -jar ${MELT_DIR}MELT.jar GroupAnalysis \
	-discoverydir $ANALYSIS_DIR \
	-w $ANALYSIS_DIR \
	-t $MEI_LIST \
	-h $GENOME_FASTA \
	-n $GENE_ANNO 
	

# Step 4 Genotype
java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
	-bamfile $WORKING_DIR$BAM1 \
	-t $MEI_LIST \
	-h $GENOME_FASTA \
	-w $ANALYSIS_DIR \
	-p $ANALYSIS_DIR 

java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
        -bamfile $WORKING_DIR$BAM2 \
        -t $MEI_LIST \
        -h $GENOME_FASTA \
        -w $ANALYSIS_DIR \
        -p $ANALYSIS_DIR 

java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
        -bamfile $WORKING_DIR$BAM3 \
        -t $MEI_LIST \
        -h $GENOME_FASTA \
        -w $ANALYSIS_DIR \
        -p $ANALYSIS_DIR 


# Make VCF
java -Xmx2G -jar ${MELT_DIR}MELT.jar MakeVCF \
	-genotypingdir $ANALYSIS_DIR \
	-h $GENOME_FASTA \
	-t $MEI_LIST \
	-w $ANALYSIS_DIR \
	-p $ANALYSIS_DIR \


