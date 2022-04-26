# Activate AWS
source /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda activate AWS

# Run this first
# aws configure 


# Run the copy
## Dry run
aws s3 sync --dryrun s3://thousandgenometriovariants.store.ubc-hpc.cloud /mnt/common/OPEN_DATA/1kG_Trio/ --include "*" --exclude ""

## Actual command

aws s3 sync --quiet s3://thousandgenometriovariants.store.ubc-hpc.cloud /mnt/common/OPEN_DATA/1kG_Trio/ --include "*" --exclude ""


