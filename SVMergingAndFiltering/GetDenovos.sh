# Set files up
SVFILE=InHouseDB_SV_50_20210827_noTRA-annotsv-candidateGenes.tsv
SVFILE_NOSPLIT=InHouseDB_SV_50_20210827_noTRA-annotsv-candidateGenes_NoSplit.tsv
SVFILE_NOSPLIT_DENOVO=InHouseDB_SV_50_20210827_noTRA-annotsv-candidateGenes_NoSplit_Denovo.tsv

# Step 1 - cut out the split annotations
grep -v -w split $SVFILE > $SVFILE_NOSPLIT

# Step 2 - Get denovosby filtering
python /mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/GetDeNovos.py \
	-I $SVFILE_NOSPLIT \
	-O $SVFILE_NOSPLIT_DENOVO


# Set files up
SVFILE=InHouseDB_MEI_15_20210827-annotsv-candidateGenes.tsv
SVFILE_NOSPLIT=InHouseDB_MEI_15_20210827-annotsv-candidateGenes_NoSplit.tsv
SVFILE_NOSPLIT_DENOVO=InHouseDB_MEI_15_20210827-annotsv-candidateGenes_NoSplit_DeNovo.tsv

# Step 1 - cut out the split annotations
grep -v -w split $SVFILE > $SVFILE_NOSPLIT

# Step 2 - Get denovosby filtering
python /mnt/scratch/Precision/EPGEN/PROCESS/MergedSV/GetDeNovos.py \
        -I $SVFILE_NOSPLIT \
        -O $SVFILE_NOSPLIT_DENOVO




