# AnnotateVariants
This is a pipeline for variant annotation in the diagnosis of rare genetic disorders. It relies on open source data and has instructions for software installs.

## Overview
1. Set-up 
	- Prepare Datasets and databases
	- Install Necessary Tools

2. Run Test

3. Run Sample




## Set-up
*THERE IS A LOT OF SETTING UP TO DO!*  
But, once you get set up, then things run nice and smooth.

1. Install Necessary software
+ bgzip and tabix
+ vt
+ vcftools/bcftools
+ snpEff
+ vcfanno
+ vcf2db
+ gemini

2. Prepare Third-party Datasets/databases  
+ All the gemini databases
+ Polyphen2  https://github.com/quinlan-lab/pathoscore/blob/master/score-sets/GRCh37/polyphen2/make.sh
+ CADD  http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz 
+ ReMM  http://remm.visze.de/files/ReMM.v0.3.1.tsv.gz
+ gnomAD 

3. Prepare in-house variant database
+ For instructions see: https://github.com/Phillip-a-richmond/BuildInHouseDB

1) PolyPhen2
GetPolyPhen2.sh is a script I got from Brent Pedersen to set up my Polyphen.

https://github.com/quinlan-lab/pathoscore/blob/master/score-sets/GRCh37/polyphen2/make.sh

The file I use inside of the polyphen2 directory is the polyphen2.txt.gz
whole_genome_SNVs.tsv.gz is the entire genome CADD score, and comes from:

3) Genomiser ReMM score 
ReMM.v0.3.1.tsv.gz
Downloaded from http://remm.visze.de/files/ReMM.v0.3.1.tsv.gz





## Run a test


