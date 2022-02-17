# Set annotSV output dir where these files are
ANNOTSV_OUTPUT_DIR=/mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/20211130_AnnotSV/

cd $ANNOTSV_OUTPUT_DIR

###########
# All SVs #
###########
##########
# Smoove #
##########

# Set files up
SVFILE=InHouseDB_SV_50_20211130_noTRA_annotsv.tsv
SVFILE_NOSPLIT=InHouseDB_SV_50_20211130_noTRA_annotsv_noSplit.tsv
SVFILE_NOSPLIT_DENOVO=InHouseDB_SV_50_20211130_noTRA_annotsv_noSplit_DeNovo.tsv

# Step 1 - cut out the split annotations
grep -v -w split $SVFILE > $SVFILE_NOSPLIT

# Step 2 - Get denovosby filtering
python /mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/GetDeNovos.py \
	-I $SVFILE_NOSPLIT \
	-O $SVFILE_NOSPLIT_DENOVO

#######################
# Candidate Genes SVs #
#######################

# SVs over candidate genes
# Set files up
SVFILE=InHouseDB_SV_50_20211130_noTRA_annotsv-candidateGenes.tsv
SVFILE_NOSPLIT=InHouseDB_SV_50_20211130_noTRA_annotsv-candidateGenes_noSplit.tsv
SVFILE_NOSPLIT_DENOVO=InHouseDB_SV_50_20211130_noTRA_annotsv-candidateGenes_noSplit_DeNovo.tsv

# Step 1 - cut out the split annotations
grep -v -w split $SVFILE > $SVFILE_NOSPLIT

# Step 2 - Get denovosby filtering
python /mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/GetDeNovos.py \
	-I $SVFILE_NOSPLIT \
	-O $SVFILE_NOSPLIT_DENOVO

#########
# Manta #
#########

# Set files up
SVFILE=InHouseDB_SV_Manta_50_20211130_noTRA_annotsv.tsv
SVFILE_NOSPLIT=InHouseDB_SV_Manta_50_20211130_noTRA_annotsv_noSplit.tsv
SVFILE_NOSPLIT_DENOVO=InHouseDB_SV_Manta_50_20211130_noTRA_annotsv_noSplit_DeNovo.tsv

# Step 1 - cut out the split annotations
grep -v -w split $SVFILE > $SVFILE_NOSPLIT

# Step 2 - Get denovosby filtering
python /mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/GetDeNovos.py \
	-I $SVFILE_NOSPLIT \
	-O $SVFILE_NOSPLIT_DENOVO

#######################
# Candidate Genes SVs #
#######################

# SVs over candidate genes
# Set files up
SVFILE=InHouseDB_SV_Manta_50_20211130_noTRA_annotsv-candidateGenes.tsv
SVFILE_NOSPLIT=InHouseDB_SV_Manta_50_20211130_noTRA_annotsv-candidateGenes_noSplit.tsv
SVFILE_NOSPLIT_DENOVO=InHouseDB_SV_Manta_50_20211130_noTRA_annotsv-candidateGenes_noSplit_DeNovo.tsv

# Step 1 - cut out the split annotations
grep -v -w split $SVFILE > $SVFILE_NOSPLIT

# Step 2 - Get denovosby filtering
python /mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/GetDeNovos.py \
	-I $SVFILE_NOSPLIT \
	-O $SVFILE_NOSPLIT_DENOVO


############
# All MEIs #
############

# MEIs over candidate genes
# Set files up
SVFILE=InHouseDB_MEI_15_20211130_annotsv.tsv
SVFILE_NOSPLIT=InHouseDB_MEI_15_20211130_annotsv_NoSplit.tsv
SVFILE_NOSPLIT_DENOVO=InHouseDB_MEI_15_20211130_annotsv_NoSplit_DeNovo.tsv

# Step 1 - cut out the split annotations
grep -v -w split $SVFILE > $SVFILE_NOSPLIT

# Step 2 - Get denovosby filtering
python /mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/GetDeNovos.py \
        -I $SVFILE_NOSPLIT \
        -O $SVFILE_NOSPLIT_DENOVO

########################
# Candidate Genes MEIs #
########################
# Set files up
SVFILE=InHouseDB_MEI_15_20211130_annotsv-candidateGenes.tsv
SVFILE_NOSPLIT=InHouseDB_MEI_15_20211130_annotsv-candidateGenes_NoSplit.tsv
SVFILE_NOSPLIT_DENOVO=InHouseDB_MEI_15_20211130_annotsv-candidateGenes_NoSplit_DeNovo.tsv

# Step 1 - cut out the split annotations
grep -v -w split $SVFILE > $SVFILE_NOSPLIT

# Step 2 - Get denovosby filtering
python /mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/GetDeNovos.py \
        -I $SVFILE_NOSPLIT \
        -O $SVFILE_NOSPLIT_DENOVO





