# Author: Phillip Richmond
# Contact: phillip.a.richmond@gmail.com
# Purpose: This is a config file for the Bam2Gemini.py pipeline script. 
# 	   It will include locations for different tools, and databases.
#	   This is an attempt at making a more generalizable pipeline system.



# Tool executables
JAVA=/opt/tools/jdk1.7.0_79/bin/java
GATKJAR=/opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar
SNPEFFJAR=/opt/tools/snpEff/snpEff.jar
BGZIP=/opt/tools/tabix/bgzip
TABIX=/opt/tools/tabix/tabix
VT=/opt/tools/vt/vt
VCFANNO=/opt/tools/vcfanno/vcfanno
BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools
VCF2DB=/opt/tools/vcf2db/vcf2db.py
GEMINI=/opt/tools/gemini/bin/gemini

# Tool configs/relevant directories
VCFANNOTOML=$ANNOTATEVARIANTS_DIR/VCFAnno/VCFANNO_Config_PlusGNOMAD_PlusInHouse_SplitByPop_gnomAD_Exome_VNX.toml
VCFANNOLUA=$ANNOTATEVARIANTS_DIR/VCFAnno/custom.lua

# Database locations
DATADIR=/mnt/causes-vnx1/DATABASES/








