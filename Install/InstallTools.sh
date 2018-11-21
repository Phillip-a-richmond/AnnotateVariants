# Start Here



## Install Conda
https://conda.io/docs/user-guide/install/linux.html

## Install Brew
http://linuxbrew.sh/

## Now get the tools you need
conda install -c bioconda vcftools
conda install -c bioconda htslib
conda install -c bioconda vcfanno 
conda install -c bioconda vt

## Other tools you need:

### VCFAnno
https://github.com/brentp/vcfanno

### VCF2DB
https://github.com/quinlan-lab/vcf2db

### GEMINI 
http://gemini.readthedocs.io/en/latest/content/installation.html

LOCATION1="/somewhere/to/store/lots/of/databases/"
LOCATION2="/somewhere/to/store/your/executable/"
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
python gemini_install.py $LOCATION1 $LOCATION2
export PATH=$PATH:$LOCATION2






