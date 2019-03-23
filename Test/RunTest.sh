# Generate the Bam2Gemini script
python ../PipelineScripts/Bam2Gemini_Development.py \
	-d /mnt/causes-vnx1/PIPELINES/AnnotateVariants/Test/ \
	-p 16 -m 40G \
	-P NA12878_Trio.ped \
	-F NA12878_Trio \
	-v VCF \
	-V NA12878_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf,NA12891_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf,NA12892_BWAmem_dupremoved_realigned_HaplotypeCaller_chr20.vcf \
	-G hg19 \
	-T Exome \
	-A /mnt/causes-vnx1/PIPELINES/AnnotateVariants/

# Run it from the command line
sh /mnt/causes-vnx1/PIPELINES/AnnotateVariants/Test/NA12878_Trio_Bam2Gemini.sh




