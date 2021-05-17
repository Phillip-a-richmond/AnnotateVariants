# Phillip Richmond
# March 7th, 2019


# This script will acquire all necessary databases, and make any changes necessary to the downloaded files

# First it will define a database directory
# Then, it will get everything that comes with GEMINI, adding it into that database directory
# Then it will acquire additional variant databases
# Then download in silico metrics at the variant level
# Then get regional annotation databases
# Then get databases on the gene level
# Lastly, it will have comments on databases which are not trivial or available for open download from online


###############################################

#####################
# Variant Databases #
#####################

#######

# gnomAD
# Genome
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz.tbi

# rename to .gz
mv gnomad.genomes.r2.1.sites.vcf.bgz gnomad.genomes.r2.1.sites.vcf.gz
mv gnomad.genomes.r2.1.sites.vcf.bgz.tbi gnomad.genomes.r2.1.sites.vcf.gz.tbi

# Exome
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz.tbi

mv gnomad.exomes.r2.1.sites.vcf.bgz gnomad.exomes.r2.1.sites.vcf.gz
mv gnomad.exomes.r2.1.sites.vcf.bgz.tbi gnomad.exomes.r2.1.sites.vcf.gz.tbi

# Final Files
# gnomad.genomes.r2.1.sites.vcf.gz
# gnomad.genomes.r2.1.sites.vcf.gz.tbi
# gnomad.exomes.r2.1.sites.vcf.gz
# gnomad.exomes.r2.1.sites.vcf.gz.tbi 

#######

# ClinVar
wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

# Index
tabix clinvar.vcf.gz

# Final Files
# clinvar.vcf.gz
# clinvar.vcf.gz.tbi


#######

# dbSNP
wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi

# Final files
# All_20180423.vcf.gz
# All_20180423.vcf.gz.tbi




############################################


#######################
# In Silico Databases #
#######################


#######

# CADD
wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz.tbi
wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/gnomad.genomes.r2.0.1.sites.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/gnomad.genomes.r2.0.1.sites.tsv.gz.tbi

# Final Files
# whole_genome_SNVs.tsv.gz
# whole_genome_SNVs.tsv.gz.tbi
# gnomad.genomes.r2.0.1.sites.tsv.gz
# gnomad.genomes.r2.0.1.sites.tsv.gz.tbi

#######

# REVEL 
wget https://rothsj06.u.hpc.mssm.edu/revel/revel_all_chromosomes.csv.zip

# unzip
unzip revel_all_chromosomes.csv.zip

# change to tsv
sed -e 's/,/\t/g' -e 's/^chr/#chr/g' revel_all_chromosomes.csv > revel_all_chromosomes.tsv

bgzip revel_all_chromosomes.tsv
tabix -b2 -s1 -e2 -c '#' revel_all_chromosomes.tsv.gz

# Final Files
# revel_all_chromosomes.tsv.gz
# revel_all_chromosomes.tsv.gz.tbi

######

# Polyphen2
# Code borrowed from B. Petersen

wget ftp://genetics.bwh.harvard.edu/pph2/whess/polyphen-2.2.2-whess-2011_12.tab.tar.bz2
mkdir -p polyphen2
cd polyphen2
tar xjvf ../polyphen-2.2.2-whess-2011_12.tab.tar.bz2
set -euo pipefail
(echo -e "#chrom\tpos\tref\talt\thdiv\thvar";

for f in polyphen-2.2.2-whess-2011_12/*features.tab; do
    set -e
    s=$(dirname $f)/$(basename $f .features.tab).scores.tab
    paste $f $s | cut -f 1,2,58,64 
    done | grep -v ^# \
    | perl -pe 's/:|\//\t/g' \
    | sed 's/^chr//; s/ //g' \
    | awk '$4 != ""' \
    | sort -k1,1 -k2,2n) | bgzip -c > polyphen2.txt.gz

tabix -b2 -e2 polyphen2.txt.gz

# Final files:
# polyphen2.txt.gz
# polyphen2.txt.gz.tbi



###################################

######################
# Regional Databases #
######################

# Conserved coding regions (CCR)

wget https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz
wget https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz

# Concatenate them
zgrep -v "^#" ccrs.xchrom.v2.20180420.bed.gz > ccrs.xchrom.v2.20180420_noheader.bed.gz
cat ccrs.autosomes.v2.20180420.bed.gz ccrs.xchrom.v2.20180420_noheader.bed.gz > ccrs.combined.v2.20180420.bed.gz

# index
tabix ccrs.combined.v2.20180420.gz

# Final files:
# ccrs.combined.v2.20180420.gz
# ccrs.combined.v2.20180420.gz.tbi


#######

# RLCRs 
wget http://tcag.ca/documents/projects/RLCRs_no_Repeat_Masker.zip

# Unzip and rename
unzip RLCRs_no_Repeat_Masker.zip
mv RLCRs_no_Repeat_Masker.txt RLCRs_no_Repeat_Masker.bed

# bgzip and index
bgzip RLCRs_no_Repeat_Masker.bed
tabix RLCRs_no_Repeat_Masker.bed.gz

# Final Files
# RLCRs_no_Repeat_Masker.bed.gz
# RLCRs_no_Repeat_Masker.bed.gz.tbi


######

# Confident Regions Platinum Genomes
wget ftp://platgene_ro@ussd-ftp.illumina.com/2017-1.0/hg19/small_variants/ConfidentRegions.bed.gz

# index
tabix ConfidentRegions.bed.gz

# Final Files
# ConfidentRegions.bed.gz
# ConfidentRegions.bed.gz.tbi



#########################################################

##########################
# Gene Based Annotations #
##########################


# gnomAD Gene Constraint
wget https://storage.googleapis.com/gnomad-public/release/2.1/ht/constraint/constraint.txt.bgz

# Unzip
mv constraint.txt.bgz constraint.txt.gz
gunzip -c constraint.txt.gz > contraint.txt

# Final files:
# contraint.txt

#######

# RVIS
wget http://genic-intolerance.org/data/GenicIntolerance_v3_12Mar16.txt

# Cut needed columns
cut -f1,2,3 GenicIntolerance_v3_12Mar16.txt > RVIS_March2016.txt

# Final Files:
# RVIS_March2016.txt


######

# HPO
wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt
wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt

# Final files
ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt
ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt

####################################################


###############################

# Non distributable Databases #

###############################

# OMIM




##############


# SpliceAI

# User must make an Illumina baseSpace account. Download this file:

# whole_genome_filtered_spliceai_scores.vcf.gz

# from this location:

# https://basespace.illumina.com/analyses/133298165/files/172857685?projectId=66029966

# Then index it

# tabix whole_genome_filtered_spliceai_scores.vcf.gz




##############

# In-house







































