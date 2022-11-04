# Get Third Party Databases
> This directory details the acquisition of third party datasets/databases for the annotation pipeline to work. Without these, there really is no pipeline.

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


##### H) Prepare in-house variant database
+ For instructions see: https://github.com/Phillip-a-richmond/BuildInHouseDB


