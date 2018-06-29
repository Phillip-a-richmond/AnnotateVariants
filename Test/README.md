# Test run for Annotating Variants

## Data
The data used here comes from the platinum trio data, for the NA12878 trio.  

To get these variant files, I ran a standard variant calling pipeline using BWAmem and GATK HaplotypeCaller.

I then Extracted only the chr20 regions, ran bgzip, and generated an index with igvtools for the chr20 variant files. These files can remain here for testing purposes.

## Analysis
The script Bam2Gemini.py is currently hard coded for the system which I work on. Eventually, I will make it more generalizable. For now, it runs with a few different input options, depending what your starting point is.

If you are starting with GVCF for multiple files, then it will run GATK GenotypeGVCFs, before proceeding. If you start with a BAM file, it will make the GVCFs, then run GenotypeGVCFs, before proceeding. Lastly, if you are starting with multiple VCFs, which isn't recommended, then it will run CombineVariants on your VCFs to make a merged VCF before proceeding.

Running Bam2Gemini.py:
```
python Bam2Gemini.py \
	-d /mnt/causes-data01/data/RICHMOND/AnnotateVariants/Test/ \
	-p 16 -m 40G \
	-P NA12878_Trio.ped \
	-F NA12878_Trio \
	-v VCF \
	-V NA12878_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf,NA12891_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf,NA12892_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf \
	-G hg19
```
If you don't already have a PED file, you can make one with MakePED.py:

```
python MakePED.py \
	--proband NA12878_BWAmem,female,affected \
	--father NA12892_BWAmem,male,unaffected \
	--mother NA12891_BWAmem,female,unaffected \
	-F NA12878_Trio \
	-O NA12878_Trio.ped 
```

Running Bam2Gemini.py will generate a script for the scheduler. Again, this is currently hard coded for Moab with my own settings, but the shell commands are generalizable. 

The script it generates will be NA12878_Trio_Bam2Gemini.sh.

Run it, and it will generate your final .db file:

NA12878_Trio.db




