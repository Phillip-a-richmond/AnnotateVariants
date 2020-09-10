# gnomAD Fetch script
# NOTE: AT the download of this database the current version is 154.
# September 4, 2020

# RefSNP VCF files for GRC (Genome Reference Consortium) human assembly
# 37 (GCF_000001405.25) and 38 (GCF_000001405.38).

# GRCH37
WORKING_DIR=/mnt/common/DATABASES/REFERENCES/GRCh37/DBSNP/
mkdir -p $WORKING_DIR
cd $WORKING_DIR

wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi

mv GCF_000001405.25.gz dbSNP_All_v154.vcf.gz
mv GCF_000001405.25.gz.tbi dbSNP_All_v154.vcf.gz.tbi

# GRCh38
WORKING_DIR=/mnt/common/DATABASES/REFERENCES/GRCh38/DBSNP/
mkdir -p $WORKING_DIR
cd $WORKING_DIR

wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz.tbi

mv GCF_000001405.38.gz dbSNP_All_v154.vcf.gz
mv GCF_000001405.38.gz.tbi dbSNP_All_v154.vcf.gz.tbi

