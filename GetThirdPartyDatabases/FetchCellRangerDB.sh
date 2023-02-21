# mm10
DatabaseDir=/mnt/common/DATABASES/
GenomeDir=$DatabaseDir/REFERENCES/mm10/
mkdir -p $GenomeDir
cd $GenomeDir

mkdir CellRanger
cd CellRanger
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
tar -zxvf refdata-gex-mm10-2020-A.tar.gz

# mm10
DatabaseDir=/mnt/common/DATABASES/
GenomeDir=$DatabaseDir/REFERENCES/GRCh38/
cd $GenomeDir

mkdir CellRanger
cd CellRanger
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz

