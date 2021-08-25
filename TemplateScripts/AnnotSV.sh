# Load BCFTools/BEDTools
source /mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/opt/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/opt/AnnotateVariantsEnvironment/

# BGZip the sub VCFs
bgzip ALU.final_comp.vcf
tabix ALU.final_comp.vcf.gz

bgzip HERVK.final_comp.vcf
tabix HERVK.final_comp.vcf.gz

bgzip SVA.final_comp.vcf
tabix SVA.final_comp.vcf.gz

bgzip LINE1.final_comp.vcf
tabix LINE1.final_comp.vcf.gz

# Concatenate the MELT files
bcftools concat -a \
       ALU.final_comp.vcf.gz \
       HERVK.final_comp.vcf.gz \
       SVA.final_comp.vcf.gz \
       LINE1.final_comp.vcf.gz \
	-o MergedMEI.vcf

export ANNOTSV=/mnt/common/Precision/AnnotSV/
$ANNOTSV/bin/AnnotSV -SVinputFile MergedMEI.vcf \
	-genomeBuild GRCh38 \
	-overlap 80 \
	-reciprocal yes \
       	-outputFile MergedMEI.annotated.tsv 

