This README details both the data generation and acquisition for databases within this directory (/project/projects/def-wyeth/DATABASES/)

Contributors:
Phil Richmond (phillip.a.richmond@gmail.com)
Robin van der Lee (rvdlee@cmmt.ubc.ca)

1) GeneHancer
(This was done on a different server, and then the files were copied here)


I downloaded the GeneHancer data dump from here:
https://genecards.weizmann.ac.il/geneloc_prev/index.shtml
Selected version 4.5 for NCBI Build 37 (nope, it's 38)
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


NOTE ROBIN: I believe the GeneHancer xlsx file has 1-based coordinates, so the start coordinates of the resulting bed files would have to be adjusted (subtract 1).



I needed to fix these:
grep -v "chrUn" GeneHancer_hg19.bed  | grep -v "_random" | sed 's/^chr//' | sortBed -i - > GeneHancer_hg19_fixed.bed



