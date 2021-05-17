#!/bin/bash
#PBS -N TIDEX915_new_PrimaryPipeline
#PBS -V
#PBS -o /mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/T246/TIDEX915new.o
#PBS -e /mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/T246/TIDEX915new.e
#PBS -m bea
#PBS -M bmodi@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=40gb
## Set the max walltime for the job
#PBS -l walltime=200:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=12
NSLOTS=$PBS_NUM_PPN
umask 0002
source /opt/tools/hpcenv.sh


SAMPLE_ID='TIDEX915_BWAmem'
GENOME_FASTA='/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa'
BWA_INDEX='/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa'
CHROM='/mnt/causes-vnx1/GENOMES/GSC/SplitByChrom/'
GENOMEFILE=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.genome

WORKING_DIR='/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/T246/'
FASTQR1='/mnt/causes-vnx2/TIDE/RAW/EXOME_TIDEX/T246/TIDEX915_1.fastq.gz'
FASTQR2='/mnt/causes-vnx2/TIDE/RAW/EXOME_TIDEX/T246/TIDEX915_2.fastq.gz'

 echo "Primary Analysis Started"
date

#FastQC
/opt/tools/FastQC/fastqc --extract $FASTQR1 $FASTQR2 -o $WORKING_DIR

#Map with BWA
/opt/tools/bwa-0.7.12/bwa mem $BWA_INDEX -t $NSLOTS -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" -M $FASTQR1 $FASTQR2 > $WORKING_DIR$SAMPLE_ID'.sam'

#Convert to binary, sort, and index
/opt/tools/samtools-1.2/samtools view -@ $NSLOTS -u -bS $WORKING_DIR$SAMPLE_ID'.sam' | /opt/tools/samtools-1.2/samtools sort -@ $NSLOTS -m 3G - $WORKING_DIR$SAMPLE_ID'.sorted'
/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

#Remove Duplicates
TMPDIR=$WORKING_DIR'picardtmp/'
mkdir $TMPDIR
/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/picard-tools-1.139/picard.jar MarkDuplicates I=$WORKING_DIR$SAMPLE_ID'.sorted.bam' O=$WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' REMOVE_DUPLICATES=false TMP_DIR=$TMPDIR M=$WORKING_DIR$SAMPLE_ID'_DuplicateResults.txt'
/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam'

#Realign
/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $NSLOTS -R $GENOME_FASTA -minReads 5 -I $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' -o $WORKING_DIR$SAMPLE_ID'_indelsites.intervals' 
/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T IndelRealigner -model USE_READS -R $GENOME_FASTA -targetIntervals $WORKING_DIR$SAMPLE_ID'_indelsites.intervals' -I $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' -o $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam'
/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam'

#Clean Up
ValidateOutput=`cat $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam_PicardSAMValidate.txt'`
rm $WORKING_DIR$SAMPLE_ID'.sam'
rm $WORKING_DIR$SAMPLE_ID'.sorted.bam'
rm $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam'
zip $WORKING_DIR$SAMPLE_ID'_R1_chastitypassed_fastqc.zip' $WORKING_DIR$SAMPLE_ID'_R1_chastitypassed_fastqc'
rm $WORKING_DIR$SAMPLE_ID'_R1_chastitypassed_fastqc'
zip $WORKING_DIR$SAMPLE_ID'_R2_chastitypassed_fastqc.zip' $WORKING_DIR$SAMPLE_ID'_R2_chastitypassed_fastqc'
rm $WORKING_DIR$SAMPLE_ID'_R2_chastitypassed_fastqc'

#SNP Calling HaplotypeCaller
/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T HaplotypeCaller \
	-nct 4 --never_trim_vcf_format_field \
	--genotyping_mode DISCOVERY \
	--standard_min_confidence_threshold_for_calling 10 \
 --standard_min_confidence_threshold_for_emitting 10 \
	--min_mapping_quality_score 0 \
	--min_base_quality_score 10 \
	--minReadsPerAlignmentStart 5 \
 --minPruning 2 \
	--pcr_indel_model NONE \
	--dbsnp /opt/tools/GATK-3.5-0-g36282e4/resources/dbsnp_138.b37.excluding_sites_after_129.vcf \
	-R $GENOME_FASTA \
	-I $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \
	-o $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf

# Fix the header for this VCF, and pipe into a bgzip then Tabix
sed 's/AD,Number=./AD,Number=R/g' $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf | \
 /opt/tools/htslib-1.5/bin/bgzip > $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf.gz 
/opt/tools/htslib-1.5/bin/tabix $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf.gz 

#NormBcftools 
/opt/tools/bcftools-1.5/bin/bcftools norm \
 --fasta-ref $GENOME_FASTA \
 --multiallelics -both \
 $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf.gz \
 --output-type z \
 $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.vcf.gz 

/opt/tools/htslib-1.5/bin/tabix $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.vcf.gz 


#Filter
/opt/tools/bcftools-1.5/bin/bcftools filter \
--include 'FORMAT/AD[1]>=3 && INFO/MQ>=20 && %QUAL>=20' \
$WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.vcf.gz \
--output-type z \
--output $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.filter.vcf.gz 

#SNP Calling HaplotypeCaller GVCFmode
/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar \
 -T HaplotypeCaller -nct 4 --emitRefConfidence GVCF \
 --standard_min_confidence_threshold_for_calling 10 \
 --standard_min_confidence_threshold_for_emitting 10 \
 --min_mapping_quality_score 0 \
 --min_base_quality_score 10 \
 --minReadsPerAlignmentStart 5 \
 --minPruning 2 \
 --pcr_indel_model NONE \
 --dbsnp /opt/tools/GATK-3.5-0-g36282e4/resources/dbsnp_138.b37.excluding_sites_after_129.vcf \
 -R $GENOME_FASTA -I $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \
 -o $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.g.vcf

#Picard Validate BAM
/opt/tools/jdk1.7.0_79/bin/java -jar  /opt/tools/picard-tools-1.139/picard.jar ValidateSamFile M=VERBOSE MO=10000000 I=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' O=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam_PicardSAMValidate.txt'

#MosDepth 

/opt/tools/mosdepth-0.2.2/mosdepth -n -q 0:1:10:20:40:60:100: -t $NSLOTS \
$WORKING_DIR$SAMPLE_ID $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam 

#Plot output 

python /mnt/causes-vnx1/RICHMOND/mosdepth/scripts/plot-dist.py \
$WORKING_DIR${SAMPLE_ID}.mosdepth.global.dist.txt \
-o $WORKING_DIR${SAMPLE_ID}.mostdepthCoverage.html \
> $WORKING_DIR${SAMPLE_ID}.mostdepthCoverage.summary.txt 


 echo "Primary Analysis Finished"
date
