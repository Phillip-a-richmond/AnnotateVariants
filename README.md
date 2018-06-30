# AnnotateVariants
#### Phillip Richmond (Phillip.A.Richmond@gmail.com)

> This is a pipeline for variant annotation in the diagnosis of rare genetic disorders. It relies on open source data and has instructions for software installs.

## Overview
1. [Pipeline Summary & Diagram](#pipeline-summary-and-diagram)

2. [Set-up](#set-Up) 
	- Prepare Datasets and databases
	- Install Necessary Tools

3. [Run Test](#run-test)

4. [Run Sample](#run-sample)


### Pipeline Summary And Diagram
This pipeline was designed by Phillip Richmond in order to analyze & prioritize variants in rare genetic disease cases. Currently, the pipeline uses the following list of software in order to accomplish this task, much to the thanks of tools produced and maintained by the lab of Aaron Quinlan:
+ GEMINI
+ VCFAnno
+ VCF2DB

Furthermore, this pipeline utilizes open source datasets within it's annotation framework, including:
+ CADD 
+ gnomAD
+ OMIM*
+ ClinVar
+ UCSC RefGene
+ Entrez Gene Summary
+ HPO Term Mapping
+ MeSHOPs
+ pLI
+ RVIS
+ FATHMM-XF
+ Eigen
+ FunSeq2
+ Platinum Genomes ConfidentRegions
+ UCSC Segmental Duplications

* OMIM requires a license for use of the API/downloadable databases, which must be applied for through their website.

Currently, the pipeline is hard coded for a specific cluster that uses the Torque-Moab scheduler. However, I will expand upon this to include other schedulers such as SLURM. Also, generalizing for software install locations, or developing a single install-script via bioconda will also be performed later in 2018.


![](https://github.com/Phillip-a-richmond/AnnotateVariants/blob/master/Figure3-NewInterpretationPipeline.png)



### Set-up
*THERE IS A LOT OF SETTING UP TO DO!*  
But, once you get set up, then things run nice and smooth.

1. Install Necessary software, details in InstallTools.sh
+ bgzip and tabix
+ vt
+ vcftools/bcftools
+ snpEff
+ vcfanno
+ vcf2db
+ gemini
+ In-house Scripts 


2. Prepare Third-party Datasets/databases, There are scripts to help do this 
+ All the gemini databases
+ Polyphen2  https://github.com/quinlan-lab/pathoscore/blob/master/score-sets/GRCh37/polyphen2/make.sh
+ CADD  http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz 
+ ReMM  http://remm.visze.de/files/ReMM.v0.3.1.tsv.gz
+ gnomAD http://gnomad.broadinstitute.org/downloads


##### A) PolyPhen2
GetPolyPhen2.sh is a script I got from Brent Pedersen to set up my Polyphen.

https://github.com/quinlan-lab/pathoscore/blob/master/score-sets/GRCh37/polyphen2/make.sh

The file I use inside of the polyphen2 directory is the polyphen2.txt.gz
whole_genome_SNVs.tsv.gz is the entire genome CADD score, and comes from:

##### B) CADD
whole_genome_SNVs.tsv.gz is the entire genome CADD score, and comes from:
http://krishna.gs.washington.edu/download/CADD/v1.3/


##### C) Genomiser ReMM score 
ReMM.v0.3.1.tsv.gz
Downloaded from http://remm.visze.de/files/ReMM.v0.3.1.tsv.gz

##### D) FATHMM
The FATHMM file comes from the annovar download (/mnt/causes-data01/data/Databases/annovar/humandb/hg19_fathmm.txt)
Which I then bgzipped and copied here.
bgzip -s1 -b2 hg19_fathmm.txt > hg19_fathmm_fromAnnovar.txt.gz

NOTE: This was giving me trouble, so I instead went directly to FATHMM to get the precomputed scores. For FATHMM-MKL:
https://github.com/HAShihab/fathmm-MKL
wget http://fathmm.biocompute.org.uk/database/fathmm-MKL_Current_zerobased.tab.gz
Then created tabix index using:
tabix -p bed fathmm-MKL_Current_zerobased.tab.gz 

For FATHMM-XF:
http://fathmm.biocompute.org.uk/fathmm-xf/#download
wget http://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_noncoding.vcf.gz
wget http://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_noncoding.vcf.gz.tbi


##### E) Eigen
wget --recursive --no-parent https://xioniti01.u.hpc.mssm.edu/v1.1/
Also, downloaded the training datasets.
https://xioniti01.u.hpc.mssm.edu/TrainingTestingDatasets/TrainingDatasets_Noncoding.zip
All found on the website here:
http://www.columbia.edu/~ii2135/download.html
python BuildSubsetShell.py -T Eigen -B xioniti01.u.hpc.mssm.edu/v1.1/VISTA_positive_sorted.bed -O ExtractEigen.sh -P Eigen_Enhancer
sh ExtractEigen.sh
Then moved them all to Eigen_VISTA_Enhancers

##### F) GeneHancer
I downloaded the GeneHancer data dump from here:
https://genecards.weizmann.ac.il/geneloc_prev/index.shtml
Selected version 4.5 for NCBI Build 37
Downloaded as xlsx onto laptop, sorted by chrom then start position
then opened VIM and copied data into it, a file named: GeneHancer.tsv
Then I cut a few columns to make the file GeneHancer.bed (cut -f1,4,5,6,9 GeneHancer.tsv > GeneHancer.bed)

I'm pretty sure (90%) that these coordinates are not in hg19 (those bastards!) 
So to resolve this I'll liftOver the bed file
For that, I needed liftOver (which I downloaded from UCSC, turns out you need a license which is free for academia and long story short...I got it)
Then I downloaded the hg38ToHg19.over.chain.gz file to lift from 38 to 19

So I renamed the GeneHancer.bed to be GeneHancer_hg38.bed
mv GeneHancer.bed GeneHancer_hg38.bed
Got rid of the header line
tail -n+2 GeneHancer_hg38.bed > GeneHancer_hg38_nohead.bed

Then I lift over using liftOver
./liftOver GeneHancer_hg38_nohead.bed hg38ToHg19.over.chain.gz GeneHancer_hg19.bed UnmappedGeneHancer.bed

In this process, it looks like I lost 666 enhancers somehow...?
Anyways: onward and upward
Then I sort using sortBed
sortBed -i GeneHancer_hg19.bed > GeneHancer_hg19.sorted.bed

And then bgzip and tabix index

##### G) LINSIGHT
Downloaded the LINSIGHT precomputed .bw file from here:
https://github.com/CshlSiepelLab/LINSIGHT
wget http://compgen.cshl.edu/%7Eyihuang/tracks/LINSIGHT.bw
Then I converted from .bw to .bedgraph using this tool from UCSC:
https://genome.ucsc.edu/goldenpath/help/bigWig.html
Found here: 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph


H) Prepare in-house variant database
+ For instructions see: https://github.com/Phillip-a-richmond/BuildInHouseDB



