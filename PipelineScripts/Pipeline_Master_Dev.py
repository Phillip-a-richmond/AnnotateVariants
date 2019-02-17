import sys, os, argparse

#####################################
# Author:   Phillip Richmond        #
# Contact:  prichmond@cmmt.ubc.ca   #
# Modified: David Arenillas (DJA)   #
#           dave@cmmt.ubc.ca        #
# Open source GNU licensing         #
#####################################

#################
#Pipeline Master#
#################

#This is the revised pipeline for analysis of Whole Genome Sequencing data.  I will maintain two pipelines within this script: One for BWA-GATK HC, and the old pipeline: Bowtie2-Samtools0.1.19
# Something to note, is that the SAMPLE ID now has the version of the pipeline you ran, NOTE: REMOVED 20190131
# Added CNV calling on Mar 20, 2017
# Added SV calling with LUMPY+CNVnator, Pindel, and merged with MetaSV September 2017
# Added ERDS in January 2018
# Added Jill's GATK-HC pipeline March 2018
# Added GSC- genome OR hg19 genome selection in April 2018
# Added MosDepth June 2018
# Migrated to VNX August 2018
# Added MToolBox, --mtoolbox [mtoolbox_config_file], November 2018 (RvdL)
# Added various coverage and metric calculations for WES data, --metrics-exome, December 2018 (RvdL)
# Major overhaul of the pipeline on January 31st 2019 by PAR
	# Exome vs. Genome
	# Reorganized a bit, cleaned up _BWAmem from the sampleIDs since it's no longer meaningful
	# Added SLURM compatability
	# Removed "masked" option
	

##################
### Initialize ###
##################


# Get the args from the user
def GetArgs():
	if len(sys.argv) < 2:
		print "Re-run with the -h option"
		print "Typical Running Command:"
		print "python Pipeline_Master.py -v new -d /mnt/causes-data04/ -s TDX-468 -1 /mnt/causes-data04/data2/RAW/EXOME_TIDEX/Family121/TIDEX468_CAGATC_L001_R1_001.fastq.gz -2 /mnt/causes-data04/data2/RAW/EXOME_TIDEX/Family121/TIDEX468_CAGATC_L001_R1_001.fastq.gz -p 16 -m 100G"
		sys.exit()
	
	
	# Redone with Argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-v","--version",help="Choose which pipeline version: old or new\nold: Bowtie2-Samtools\nnew:BWA-GATK HC\nIf you put down 'neither', then it will not run the primary pipeline",required=True)
	parser.add_argument("-d","--workingDir",help="Working directory on your filesystem, full filepath",required=True)
	parser.add_argument("-s","--sampleID",help="What is the sample ID for these input files, e.g. TIDEX0001",required=True)
	parser.add_argument("-p","--processors",help="Choose the number of processors for this job, e.g. -p 16",type=int)
	parser.add_argument("-m","--memory",help="Choose the memory needed for this job, e.g. -m 100G",required=True)
	parser.add_argument("-1","--inputR1",help="The R1 fastq file",required=True)
	parser.add_argument("-2","--inputR2",help="The R2 fastq file",required=True)
	parser.add_argument("-n","--nodeName",help="20190131-Deprecated!! Name of the node you want to submit to, options are: cfricauseshpc02,cfricauseshpc03,cfricauseshpc04,cfricauseshpc06")
	parser.add_argument("--cnv",help="Run CNV calling and annotation",action='store_true',default=False)
	parser.add_argument("-R","--readlength",help="The length of the reads in integer form",default=100)
	parser.add_argument("-T","--Type",help="Exome || Genome",required=True)
	parser.add_argument("-S","--scheduler",help="Which scheduler you want to submit to.  This will determine the format of the shell script. Options: PBS || SGE",type=str)
	parser.add_argument("-G","--GENOME",help="Which Genome version do you want to use? Options are GSC || hg19",required=True)
	parser.add_argument("--mtoolbox",help="Provide an MToolBox config file here to perform a mitochondrial variant analysis, e.g. /path/to/AnnotateVariants/MToolBox_config_files/MToolBox_rCRS_config_with_markdup_and_indelrealign_RvdL.sh",default="")
	parser.add_argument("--metrics-exome",help="Calculate exome coverage and other metrics using Mosdepth and Picard CalculateHsMetrics, assuming the Agilent_SureSelect_Human_All_Exon_V4 capture kit was used",action='store_true',default=False)
	parser.add_argument("--sv",help="Run SV calling and annotation",action='store_true',default=False)
	parser.add_argument("-E","--Email",help="Email address",type=str)
	args = parser.parse_args()
	return args



#########################################################################################


# Now define these functions, which when called, they will write out the command to the shell
# script.  This is an easy way to control what goes into the script, so if I have to run
# only fastQC on all the files, then I can just specify that, and call that command

###### Old Pipeline, Legacy Commands ######

# Part of the old pipeline, maintaining for legacy purposes
def Bowtie2():
	#Map with Bowtie2
	shellScriptFile.write("\n#Map with Bowtie2\n")
	shellScriptFile.write("/opt/tools/bowtie2-2.2.6/bowtie2  -x $BOWTIE2_INDEX -1 $FASTQR1 -2 $FASTQR2 -S $WORKING_DIR${SAMPLE_ID}.sam  -p $NSLOTS --very-sensitive -X 1000 --met-stderr --rg-id $SAMPLE_ID --rg \"SM:$SAMPLE_ID\tPL:illumina\" 2> $WORKING_DIR${SAMPLE_ID}.stderr\n\n")
	shellScriptFile.write("\n#Fix Read Names\n")
	shellScriptFile.write("#python /mnt/causes-vnx1/PipelineControl/FixSamReadNames.py $WORKING_DIR${SAMPLE_ID}.sam $WORKING_DIR${SAMPLE_ID}fixed.sam\n")
	shellScriptFile.write("#mv $WORKING_DIR$SAMPLE_ID\'fixed.sam\' $WORKING_DIR${SAMPLE_ID}.sam\n")

#SNP calling mpileup
def VarCall(shellScriptFile):
	shellScriptFile.write("\n#SNPcalling Mpileup\n")
	shellScriptFile.write("#/opt/tools/samtools-1.2/samtools mpileup -Bgf $GENOME_FASTA $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\' | /opt/tools/bcftools-1.2/bcftools call -vc - > $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned_mpileup.bcf\'\n")
	shellScriptFile.write("#/opt/tools/bcftools-1.2/vcfutils.pl varFilter -Q 30 -a 2 $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned_mpileup.bcf\' | awk \'(match ($1, \"##\") || $6 > 30)\' > $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned_mpileup.vcf\'\n")

#old mpileup
	shellScriptFile.write("\n#SNP calling Mpileupv0.1.19\n")
	shellScriptFile.write("/opt/tools/samtools-0.1.19/samtools mpileup -Bgf $GENOME_FASTA $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' | /opt/tools/samtools-0.1.19/bcftools/bcftools view -gvc - > $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_v0.1.19mpileup.bcf'")
	shellScriptFile.write("\n/opt/tools/samtools-0.1.19/bcftools/vcfutils.pl varFilter -Q20 -a 3 $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_v0.1.19mpileup.bcf' | awk '(match ($1,\"##\") || $6 > 30)' > $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_v0.1.19mpileup.vcf'\n")

# Part of the old pipeline
def IndelRemoval(shellScriptFile):
	shellScriptFile.write("\n#Indel Removal\n")
	shellScriptFile.write("/opt/tools/tabix-0.2.6/bgzip -c $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_v0.1.19mpileup.vcf' > $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_v0.1.19mpileup.vcf.gz' \n")
	shellScriptFile.write("/opt/tools/tabix-0.2.6/tabix $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_v0.1.19mpileup.vcf.gz'\n")
	shellScriptFile.write("/opt/tools/bcftools-1.2/bcftools filter -e \"TYPE='indel' && IS=1\"  $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_v0.1.19mpileup.vcf.gz' > $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_v0.1.19mpileup_rmindelIS1.vcf' \n")

# Added DJA 2016/06/15
def SNVMetrics(shellScriptFile):
	shellScriptFile.write("\n#SNV Metrics\n")
	shellScriptFile.write("java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T VariantEval --eval $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned_v0.1.19mpileup_rmindelIS1.vcf\' --out $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned_v0.1.19mpileup_rmindelIS1.vcf_metrics\' --dbsnp /opt/tools/GATK-3.4-46/resources/dbsnp_138.hg19_noM.vcf --intervals /mnt/causes-vnx1/mwenifumbo/ExomeCapture/hg19_RefSeq_050116_CodingExons_nohap.interval_list -R $GENOME_FASTA \n")

def SummaryStats(shellScriptFile):
	#Parse Stats
	shellScriptFile.write("\n#Parse Stats\n")
	shellScriptFile.write("python /mnt/causes-vnx1/PipelineControl/ParseStats.py %s %s %s %s > $WORKING_DIR$SAMPLE_ID'_summaryStats.txt'\n"%(sampleID,workingDir,R1fastq,R2fastq))

def Platypus(shellScriptFile):
	shellScriptFile.write("\n#SNP Calling Platypus\n")
	shellScriptFile.write("python /opt/tools/Platypus-0.8.1/Platypus.py callVariants --nCPU=$NSLOTS --refFile=$GENOME_FASTA --bamFiles=$WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\' --output=$WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned_Platypus.vcf\'\n")
	shellScriptFile.write("\n#Sort platypus file\n")
	shellScriptFile.write("mkdir $WORKING_DIR'tmp/' \n")
	shellScriptFile.write("/opt/tools/vcftools-0.1.14/bin/vcf-sort -c -p $NSLOTS -t $WORKING_DIR'tmp/' $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_Platypus.vcf' > $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_Platypus_sorted.vcf'\n")
	shellScriptFile.write("mv $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_Platypus_sorted.vcf'  $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_Platypus.vcf'\n")

###### Either new pipeline, or shared pipeline commands (both old and new) #######

# FastQC to generate quality report. No interpretation of the quality report though
def FastQC(shellScriptFile):
	shellScriptFile.write("\n#FastQC\n")
	shellScriptFile.write("/opt/tools/FastQC/fastqc --extract $FASTQR1 $FASTQR2 -o $WORKING_DIR\n")

#Map with BWA
def BWA(shellScriptFile):
	shellScriptFile.write("\n#Map with BWA\n")
	shellScriptFile.write("/opt/tools/bwa-0.7.12/bwa mem $BWA_INDEX -t $NSLOTS -R \"@RG\\tID:$SAMPLE_ID\\tSM:$SAMPLE_ID\\tPL:illumina\" -M $FASTQR1 $FASTQR2 > $WORKING_DIR$SAMPLE_ID'.sam'\n\n")

#Convert to binary, sort, and index
def Sam2Bam(shellScriptFile):
	shellScriptFile.write("#Convert to binary, sort, and index\n")
	shellScriptFile.write("/opt/tools/samtools-1.2/samtools view -@ $NSLOTS -u -bS $WORKING_DIR$SAMPLE_ID\'.sam\' | /opt/tools/samtools-1.2/samtools sort -@ $NSLOTS -m 3G - $WORKING_DIR$SAMPLE_ID\'.sorted\'\n")
	shellScriptFile.write("/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID\'.sorted.bam\'\n")

# Duplicate marking
def DupRemove(shellScriptFile):
	shellScriptFile.write("\n#Remove Duplicates\n")
	shellScriptFile.write("TMPDIR=$WORKING_DIR'picardtmp/'\n")
	shellScriptFile.write("mkdir $TMPDIR\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/picard-tools-1.139/picard.jar MarkDuplicates I=$WORKING_DIR$SAMPLE_ID\'.sorted.bam\' O=$WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\' REMOVE_DUPLICATES=false TMP_DIR=$TMPDIR M=$WORKING_DIR$SAMPLE_ID\'_DuplicateResults.txt\'\n")
	shellScriptFile.write("/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\'\n")

# GATK Realignment
def Realign(shellScriptFile):
	shellScriptFile.write("\n#Realign\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $NSLOTS -R $GENOME_FASTA -minReads 5 -I $WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_indelsites.intervals\' \n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T IndelRealigner -model USE_READS -R $GENOME_FASTA -targetIntervals $WORKING_DIR$SAMPLE_ID\'_indelsites.intervals\' -I $WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\'\n")
	shellScriptFile.write("/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\'\n")

# Added DJA 2016/06/15
def ReadMetrics(shellScriptFile):
	shellScriptFile.write("\n#Read Metrics\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/picard-tools-1.139/picard.jar CollectMultipleMetrics INPUT=$WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\' OUTPUT=$WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted\' REFERENCE_SEQUENCE=$GENOME_FASTA\n")

# Functional
def ValidateSAM(shellScriptFile):
	shellScriptFile.write("\n#Picard Validate BAM\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar  /opt/tools/picard-tools-1.139/picard.jar ValidateSamFile M=VERBOSE MO=10000000 I=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' O=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam_PicardSAMValidate.txt'\n")

def CleanUp(shellScriptFile):
	#Clean Up
	shellScriptFile.write("\n#Clean Up\n")
	shellScriptFile.write("ValidateOutput=`cat $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam_PicardSAMValidate.txt'`\n")
	shellScriptFile.write("rm $WORKING_DIR$SAMPLE_ID\'.sam\'\n")
	shellScriptFile.write("rm $WORKING_DIR$SAMPLE_ID\'.sorted.bam\'\n")
	shellScriptFile.write("rm $WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\'\n")
	shellScriptFile.write("zip $WORKING_DIR$SAMPLE_ID\'_R1_chastitypassed_fastqc.zip\' $WORKING_DIR$SAMPLE_ID\'_R1_chastitypassed_fastqc\'\n")
	shellScriptFile.write("rm $WORKING_DIR$SAMPLE_ID\'_R1_chastitypassed_fastqc\'\n")
	shellScriptFile.write("zip $WORKING_DIR$SAMPLE_ID\'_R2_chastitypassed_fastqc.zip\' $WORKING_DIR$SAMPLE_ID\'_R2_chastitypassed_fastqc\'\n")
	shellScriptFile.write("rm $WORKING_DIR$SAMPLE_ID\'_R2_chastitypassed_fastqc\'\n")


#SNP calling with GATK HaplotypeCaller, Jmwenifumbo's code here. Also includes BCFtools normalization and filter
def GATK_HaplotypeCaller(shellScriptFile):
	shellScriptFile.write("\n#SNP Calling HaplotypeCaller\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T HaplotypeCaller \\\n")
	shellScriptFile.write("	-nct 4 --never_trim_vcf_format_field \\\n")
	shellScriptFile.write("	--genotyping_mode DISCOVERY \\\n")
	shellScriptFile.write("	--standard_min_confidence_threshold_for_calling 10 \\\n")
	shellScriptFile.write(" --standard_min_confidence_threshold_for_emitting 10 \\\n")
	shellScriptFile.write("	--min_mapping_quality_score 0 \\\n")
	shellScriptFile.write("	--min_base_quality_score 10 \\\n")
	shellScriptFile.write("	--minReadsPerAlignmentStart 5 \\\n")
	shellScriptFile.write(" --minPruning 2 \\\n")
	shellScriptFile.write("	--pcr_indel_model NONE \\\n")
	shellScriptFile.write("	--dbsnp /opt/tools/GATK-3.5-0-g36282e4/resources/dbsnp_138.b37.excluding_sites_after_129.vcf \\\n")
	shellScriptFile.write("	-R $GENOME_FASTA \\\n")
	shellScriptFile.write("	-I $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
	shellScriptFile.write("	-o $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf\n")
	shellScriptFile.write("\n# Fix the header for this VCF, and pipe into a bgzip then Tabix\n")
	shellScriptFile.write("sed 's/AD,Number=./AD,Number=R/g' $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf | \\\n")
	shellScriptFile.write(" /opt/tools/htslib-1.5/bin/bgzip > $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf.gz \n")
	shellScriptFile.write("/opt/tools/htslib-1.5/bin/tabix $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf.gz \n")
	shellScriptFile.write("\n#NormBcftools \n")
	shellScriptFile.write("/opt/tools/bcftools-1.5/bin/bcftools norm \\\n")
	shellScriptFile.write(" --fasta-ref $GENOME_FASTA \\\n")
	shellScriptFile.write(" --multiallelics -both \\\n")
	shellScriptFile.write(" $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf.gz \\\n")
	shellScriptFile.write(" --output-type z \\\n")
	shellScriptFile.write(" $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.vcf.gz \n\n")
	shellScriptFile.write("/opt/tools/htslib-1.5/bin/tabix $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.vcf.gz \n\n")
	shellScriptFile.write("\n#Filter\n")
	shellScriptFile.write("/opt/tools/bcftools-1.5/bin/bcftools filter \\\n")
	shellScriptFile.write("--include 'FORMAT/AD[1]>=3 && INFO/MQ>=20 && %QUAL>=20' \\\n")
	shellScriptFile.write("$WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.vcf.gz \\\n")
	shellScriptFile.write("--output-type z \\\n")
	shellScriptFile.write("--output $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.filter.vcf.gz \n")

#SNP calling, creating a GVCF which is used in multi-sample variant calling (Bam2Gemini.sh part of the pipeline)
def GATK_HaplotypeCallerGVCF(shellScriptFile):
	shellScriptFile.write("\n#SNP Calling HaplotypeCaller GVCFmode\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar \\\n")
	shellScriptFile.write(" -T HaplotypeCaller -nct 4 --emitRefConfidence GVCF \\\n")
        shellScriptFile.write(" --standard_min_confidence_threshold_for_calling 10 \\\n")
        shellScriptFile.write(" --standard_min_confidence_threshold_for_emitting 10 \\\n")
        shellScriptFile.write(" --min_mapping_quality_score 0 \\\n")
        shellScriptFile.write(" --min_base_quality_score 10 \\\n")
        shellScriptFile.write(" --minReadsPerAlignmentStart 5 \\\n")
        shellScriptFile.write(" --minPruning 2 \\\n")
        shellScriptFile.write(" --pcr_indel_model NONE \\\n")
        shellScriptFile.write(" --dbsnp /opt/tools/GATK-3.5-0-g36282e4/resources/dbsnp_138.b37.excluding_sites_after_129.vcf \\\n")
	shellScriptFile.write(" -R $GENOME_FASTA -I $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
	shellScriptFile.write(" -o $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.g.vcf\n")

# GATK Depth of Coverage Assessment
def GATK_Coverage(shellScriptFile):
	shellScriptFile.write("\n#GATK Depth Of Coverage\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $GENOME_FASTA -I $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned_Coverage'\n")


# MosDepth depth of coverage assessment
def MosDepth_WGS(shellScriptFile):
	shellScriptFile.write("\n#MosDepth \n\n")
	shellScriptFile.write("/opt/tools/mosdepth-0.2.2/mosdepth -n -q 0:1:10:20:40:60:100: -t $NSLOTS \\\n")
	shellScriptFile.write("$WORKING_DIR$SAMPLE_ID $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	shellScriptFile.write("\n#Plot output \n\n")
	shellScriptFile.write("python /opt/tools/mosdepth-0.2.2/plot-dist.py \\\n")
	shellScriptFile.write("$WORKING_DIR${SAMPLE_ID}.mosdepth.global.dist.txt \\\n")
	shellScriptFile.write("-o $WORKING_DIR${SAMPLE_ID}.mostdepthCoverage.html \\\n")
	shellScriptFile.write("> $WORKING_DIR${SAMPLE_ID}.mostdepthCoverage.summary.txt \n\n")

def MosDepth_WES(shellScriptFile):
	shellScriptFile.write("\n# Run MosDepth \n")
	shellScriptFile.write("SAMPLE=\'%s\'\n"%sampleID)
	shellScriptFile.write("EXOME_CAPTURE_BED=/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/Agilent_SureSelect_Human_All_Exon_V4/S03723314_Covered_chrnameswithoutchr.bed \n")
	shellScriptFile.write("MOSDEPTH_PATH=/opt/tools/mosdepth-0.2.2/ \n")
	shellScriptFile.write(" \n")
	shellScriptFile.write("$MOSDEPTH_PATH/mosdepth -t 4 -n -b $EXOME_CAPTURE_BED $METRICS_WORKING_DIR/$SAMPLE $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	shellScriptFile.write("python $MOSDEPTH_PATH/plot-dist.py -o $METRICS_WORKING_DIR/mosdepth_coverage.html $METRICS_WORKING_DIR/*region.dist.txt \n")
	shellScriptFile.write("COV_OUT=$METRICS_WORKING_DIR/mean_coverage_mosdepth.txt \n")
	shellScriptFile.write("python $MOSDEPTH_PATH/mean_coverage_mosdepth.py $METRICS_WORKING_DIR/*regions.bed.gz > $COV_OUT \n")
	shellScriptFile.write("head -1000 $COV_OUT \n")
	shellScriptFile.write(" \n")

def Picard_HSMETRICS(shellScriptFile):
	shellScriptFile.write("\n# Run Picard CalculateHsMetrics \n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/picard-tools-1.139/picard.jar CalculateHsMetrics \\\n")
	shellScriptFile.write("\t I=$WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
	shellScriptFile.write("\t O=$METRICS_WORKING_DIR${SAMPLE_ID}_picard_hsmetrics.txt \\\n")
	shellScriptFile.write("\t R=$GENOME_FASTA \\\n")
	shellScriptFile.write("\t BAIT_INTERVALS=$EXOME_CAPTURE_INTERVAL \\\n")
	shellScriptFile.write("\t TARGET_INTERVALS=$EXOME_CAPTURE_INTERVAL \n")
	shellScriptFile.write(" \n")


# MToolBox Mitochondrial analysis
def MToolBox(shellScriptFile):
	shellScriptFile.write("\n# Run MToolBox for mitochondrial variant analysis \n")
	shellScriptFile.write("SAMPLE=\'%s\'\n"%sampleID)
	shellScriptFile.write("MTOOLBOX_PATH=/opt/tools/MToolBox-1.0/ \n")
	shellScriptFile.write("#MTOOLBOX_PATH=/mnt/home/BCRICWH.LAN/rvanderlee/MToolBox-1.1/ \n")
	shellScriptFile.write("PATH=$MTOOLBOX_PATH/MToolBox/:$MTOOLBOX_PATH:$PATH \n")
	shellScriptFile.write("MTOOLBOX_WORKING_DIR=$WORKING_DIR/MToolBox_${SAMPLE}/ \n")
	shellScriptFile.write("mkdir -p $MTOOLBOX_WORKING_DIR/ \n")
	shellScriptFile.write(" \n")
	shellScriptFile.write("MTOOLBOX_CONFIG_FILE_ORIGINAL=\'%s\'\n"%mtoolboxConfigFile)
	shellScriptFile.write("MTOOLBOX_CONFIG_FILE_ORIGINAL_BASENAME=$(basename $MTOOLBOX_CONFIG_FILE_ORIGINAL) \n")
	shellScriptFile.write("MTOOLBOX_CONFIG_FILE=$MTOOLBOX_WORKING_DIR/$MTOOLBOX_CONFIG_FILE_ORIGINAL_BASENAME \n")
	shellScriptFile.write(" \n")
	shellScriptFile.write("#edit the MToolBox config file template so that is specifies the MToolBox results/working directory for the current analysis \n")
	shellScriptFile.write("cp $MTOOLBOX_CONFIG_FILE_ORIGINAL $MTOOLBOX_CONFIG_FILE \n")
	shellScriptFile.write("sed -i \"s#^output_name\=\.#output_name=$MTOOLBOX_WORKING_DIR#\" $MTOOLBOX_CONFIG_FILE \n" )
	shellScriptFile.write(" \n")
	shellScriptFile.write("#link the raw fastq files to the MToolBox working directory, and name them as required by MToolBox: \<sample\_name\>.R1.fastq, \<sample\_name\>.R2.fastq \n")
	shellScriptFile.write("ln -sf ${FASTQR1} $MTOOLBOX_WORKING_DIR/${SAMPLE}.R1.fastq.gz \n")
	shellScriptFile.write("ln -sf ${FASTQR2} $MTOOLBOX_WORKING_DIR/${SAMPLE}.R2.fastq.gz \n")
	shellScriptFile.write(" \n")
	shellScriptFile.write("echo \"Changing working directory to $MTOOLBOX_WORKING_DIR and running MToolBox...\" \n")
	shellScriptFile.write("PWD_CURRENT=`pwd` \n")
	shellScriptFile.write("cd $MTOOLBOX_WORKING_DIR \n")
	shellScriptFile.write("$MTOOLBOX_PATH/MToolBox/MToolBox.sh -i ${MTOOLBOX_CONFIG_FILE} \n")
	shellScriptFile.write("echo \"Changing working directory to back to $PWD_CURRENT...\" \n")
	shellScriptFile.write("cd $PWD_CURRENT \n")
	shellScriptFile.write(" \n")


# Structural variant calling

# CNV Calling with CNVnator
def CNVNATOR(shellScriptFile):
	shellScriptFile.write("#Running CNVNATOR Windowsize 100\n\n")
	shellScriptFile.write("#Defining Variables\n")
	shellScriptFile.write("BAM=${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	shellScriptFile.write("WIN=100\n")
	shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR${BAM}.root -genome $GENOME_FASTA -tree $WORKING_DIR$BAM \n")
	shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -his $WIN -d $CHROM \n")
	shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -stat $WIN \n")
	shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -partition $WIN \n")
	shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -call $WIN > $WORKING_DIR/${SAMPLE_ID}_CNVnatorCall_$WIN \n")

	shellScriptFile.write("# Preparing for use within Lumpy\n")
	shellScriptFile.write("/opt/tools/lumpy-0.2.11/scripts/cnvanator_to_bedpes.py -c $WORKING_DIR${SAMPLE_ID}_CNVnatorCall_$WIN -b 600 --del_o $WORKING_DIR${SAMPLE_ID}.del.$WIN.bedpe --dup_o $WORKING_DIR${SAMPLE_ID}.dup.$WIN.bedpe \n")
	shellScriptFile.write("/opt/tools/lumpy-0.2.11/scripts/bedpe_sort.py -b $WORKING_DIR${SAMPLE_ID}.del.$WIN.bedpe -g $GENOMEFILE > $WORKING_DIR${SAMPLE_ID}.del.$WIN.sorted.bedpe \n")
	shellScriptFile.write("/opt/tools/lumpy-0.2.11/scripts/bedpe_sort.py -b $WORKING_DIR${SAMPLE_ID}.dup.$WIN.bedpe -g $GENOMEFILE > $WORKING_DIR${SAMPLE_ID}.dup.$WIN.sorted.bedpe \n")

	shellScriptFile.write("#Running CNVNATOR Windowsize 1000\n\n")
        shellScriptFile.write("WIN=1000\n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -his $WIN -d $CHROM \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -stat $WIN \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -partition $WIN \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -call $WIN > $WORKING_DIR/${SAMPLE_ID}_CNVnatorCall_$WIN \n")

#CNV Calling with ERDS
def ERDS(shellScriptFile):
	shellScriptFile.write("\n# ERDS \n\n")
	shellScriptFile.write("perl /opt/tools/erds1.1/erds_pipeline.pl \\\n")
	shellScriptFile.write("\t-b $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
	shellScriptFile.write("\t-v $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf \\\n")
	shellScriptFile.write("\t-o $WORKING_DIR${SAMPLE_ID}_ERDS/ \\\n")
	shellScriptFile.write("\t-r $GENOME_FASTA \\\n")

# CNV calling with CANVAS, also includes an annovar call
def CANVAS(shellScriptFile):
	shellScriptFile.write("#Running CANVAS\n#Redefining Variables\n")
	shellScriptFile.write("BAMFILE=$SAMPLE_ID'_dupremoved_realigned'\n")
	shellScriptFile.write("VCFFILE=$SAMPLE_ID'_dupremoved_realigned_Platypus.vcf'\n")
	shellScriptFile.write("CANVAS_DIR=$WORKING_DIR$SAMPLE_ID'_CANVAS/'\n")
	shellScriptFile.write("CANVAS_VCFGZ='CNV.vcf.gz'\n")
	shellScriptFile.write("CANVAS_VCF=$SAMPLE_ID'_CANVAS_CNV.vcf'\n")
	shellScriptFile.write("CANVAS_VCF_FILTERED=$SAMPLE_ID'_CANVAS_CNV_filtered.vcf'\n")
	shellScriptFile.write("CANVAS_AVINPUT=$SAMPLE_ID'_CANVAS.avinput'\n")
	shellScriptFile.write("CANVAS_LOH_AVINPUT=$SAMPLE_ID'_CANVAS_LOH.avinput'\n")
	shellScriptFile.write("CANVAS_ANNOTATED=$SAMPLE_ID'_CANVAS_annovar'\n")
	shellScriptFile.write("CANVAS_BED=$SAMPLE_ID'_CANVAS_CNV.bed'\n")

	shellScriptFile.write("\n#Make the working directory for CANVAS\n")
	shellScriptFile.write("mkdir $CANVAS_DIR\n")
	shellScriptFile.write("#Run CANVAS\n")
	shellScriptFile.write("mono /opt/tools/Canvas/Canvas.exe Germline-WGS -b $WORKING_DIR$BAMFILE --b-allele-vcf=$WORKING_DIR$VCFFILE -o $CANVAS_DIR -r /mnt/causes-vnx1/GENOMES/hg19/hg19_CANVAS_kmer.fa -g /mnt/causes-vnx1/GENOMES/hg19/ -f /mnt/causes-vnx1/GENOMES/hg19/hg19_CANVAS_Filter13.bed -n $SAMPLE_ID\n\n")

	shellScriptFile.write("#Get rid of the reference calls (also unzips)\n")
	shellScriptFile.write("zcat $CANVAS_DIR$CANVAS_VCFGZ | grep -v ':REF:' > $CANVAS_DIR$CANVAS_VCF\n")
	shellScriptFile.write("#Filter using custom script\n")
	shellScriptFile.write("python /opt/tools/VariantAnnotation/FilterCanvas.py $CANVAS_DIR$CANVAS_VCF $CANVAS_DIR$CANVAS_VCF_FILTERED $CANVAS_DIR$CANVAS_AVINPUT $CANVAS_DIR$CANVAS_LOH_AVINPUT\n\n")
	shellScriptFile.write("#Copy to a bed file (nothing special)\n")
	shellScriptFile.write("cp $CANVAS_DIR$CANVAS_AVINPUT $CANVAS_DIR$CANVAS_BED\n\n")
	shellScriptFile.write("#Run annovar\n")
	shellScriptFile.write("/mnt/causes-vnx1/DATABASES/annovar/new_table_annovar.pl --otherinfo --buildver hg19 --protocol refgene,dgvMerged_commafix --operation g,r --argument '','-minqueryfrac 0.8 --colsWanted 2&3&4&10&17&18&19&21' \\\n")
	shellScriptFile.write("$CANVAS_DIR$CANVAS_AVINPUT /mnt/causes-vnx1/DATABASES/annovar/humandb -out $WORKING_DIR$CANVAS_ANNOTATED\n")

# SV Calling with Pindel
def Pindel(shellScriptFile):
	shellScriptFile.write("# Running Pindel \n\n")
	shellScriptFile.write("#Defining Variables\n")
        shellScriptFile.write("BAM=${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	shellScriptFile.write("\n##Generate Empirical insert size stats\n")
        shellScriptFile.write("/opt/tools/samtools-1.2/samtools view -r $SAMPLE_ID $WORKING_DIR$BAM | tail -n+1000000 | python /opt/tools/lumpy/pairend_distro.py -r %s -X 4 -N 1000000 -o $WORKING_DIR${BAM}.histo > $WORKING_DIR${BAM}.insertStats\n"%args.readlength)
        shellScriptFile.write("MEAN=`cat $WORKING_DIR${BAM}.insertStats | sed -E 's/\s+/,/' | cut -d, -f1 | sed -E 's/mean://' | xargs printf \"%.0f\"`\n")
	shellScriptFile.write("python /mnt/causes-vnx1/PipelineControl/processingpipeline/pindel_config.py $WORKING_DIR$BAM $MEAN $WORKING_DIR$SAMPLE_ID \n")
	shellScriptFile.write("/opt/tools/pindel-0.2.5b6/pindel --number_of_threads $NSLOTS \\\n")
	shellScriptFile.write("-f $GENOME_FASTA \\\n")
	shellScriptFile.write("-i $WORKING_DIR${BAM}_config.txt \\\n")
	shellScriptFile.write("-c ALL -o $WORKING_DIR$SAMPLE_ID \\\n")
	shellScriptFile.write("-M 4 \\\n")
	shellScriptFile.write("-N -x 3 \\\n")
	shellScriptFile.write("-I false -t false -r false \\\n")
	shellScriptFile.write("-J /mnt/causes-vnx1/DATABASES/hg19_centromeres_telomeres.bed  \n\n")

# Merging SVs with MetaSV
def MetaSV(shellScriptFile):
	shellScriptFile.write("\n# MetaSV \n\n")
	shellScriptFile.write("BAM=${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	shellScriptFile.write("run_metasv.py --bam $WORKING_DIR/$BAM --reference $GENOME_FASTA \\\n")
	shellScriptFile.write("--enable_per_tools_output \\\n")
	shellScriptFile.write("--pindel_native $WORKING_DIR/${SAMPLE_ID}_D $WORKING_DIR/${SAMPLE_ID}_TD \\\n")
	shellScriptFile.write("--lumpy_vcf $WORKING_DIR${SAMPLE_ID}_Lumpy.vcf \\\n")
	shellScriptFile.write("--cnvnator_native ${SAMPLE_ID}_CNVnatorCall_1000 \\\n")
	shellScriptFile.write("--filter_gaps --sample $WORKING_DIR$SAMPLE_ID --num_threads $NSLOTS \\\n")
	shellScriptFile.write("--outdir $WORKING_DIR${SAMPLE_ID}MetaSVout \\\n")
	shellScriptFile.write("--workdir $WORKING_DIR${SAMPLE_ID}MetaSVwork \\\n")
	shellScriptFile.write("--spades /opt/tools/SPAdes-3.10.1/bin/spades.py \\\n")
	shellScriptFile.write("--age /opt/tools/AGE/age_align \\\n")

# SV Calling Deprecated
def Lumpy(shellScriptFile):
	shellScriptFile.write("\n\n#Lumpy CNV Caller\n")
	shellScriptFile.write("\n##PreProcess\n")
	shellScriptFile.write("BAM_FILE=$SAMPLE_ID'_dupremoved_realigned'\n")
	shellScriptFile.write("/opt/tools/samtools-1.2//samtools view -@ $NSLOTS -b -F 1294 $WORKING_DIR$BAM_FILE'.sorted.bam' > $WORKING_DIR$BAM_FILE'_discordants.bam'\n")
	shellScriptFile.write("/opt/tools/samtools-1.2/samtools view -@ $NSLOTS -h $WORKING_DIR$BAM_FILE'.sorted.bam' | python /opt/tools/lumpy/extractSplitReads_BwaMem.py -i stdin | /opt/tools/samtools-1.2/samtools view -@ $NSLOTS -Sb - > $WORKING_DIR$BAM_FILE'_splitters.bam'\n")
	shellScriptFile.write("\n##Sort the splitters and discordants\n")
	shellScriptFile.write("/opt/tools/samtools-1.2/samtools sort -@ $NSLOTS $WORKING_DIR$BAM_FILE'_discordants.bam' $WORKING_DIR$BAM_FILE'_discordants.sorted'\n")
	shellScriptFile.write("/opt/tools/samtools-1.2/samtools sort -@ $NSLOTS $WORKING_DIR$BAM_FILE'_splitters.bam' $WORKING_DIR$BAM_FILE'_splitters.sorted'\n")
	shellScriptFile.write("\n##Generate Empirical insert size stats\n")
	shellScriptFile.write("/opt/tools/samtools-1.2/samtools view -r $SAMPLE_ID $WORKING_DIR$BAM_FILE'.sorted.bam' | tail -n+1000000 | python /opt/tools/lumpy/pairend_distro.py -r %s -X 4 -N 1000000 -o $WORKING_DIR$BAM_FILE'.histo' > $WORKING_DIR$BAM_FILE'.insertStats'\n"%args.readlength)
	shellScriptFile.write("MEAN=`cat $WORKING_DIR$BAM_FILE'.insertStats' | sed -E 's/\s+/,/' | cut -d, -f1 | sed -E 's/mean://' | xargs printf \"%.0f\"`\n")
	shellScriptFile.write("STDEV=`cat $WORKING_DIR$BAM_FILE'.insertStats' | sed -E 's/\s+/,/' | cut -d, -f2 | sed -E 's/stdev://' | xargs printf \"%.0f\"`\n")
	shellScriptFile.write("echo \"$STDEV\"\n")
	shellScriptFile.write("echo \"$MEAN\"\n")
	shellScriptFile.write("\n##Run LUMPY\n")
	shellScriptFile.write("/opt/tools/lumpy/lumpy -mw 4 -tt 1.0 \\\n")
	shellScriptFile.write("-bedpe bedpe_file:$WORKING_DIR${SAMPLE_ID}.dup.100.sorted.bedpe,weight:3,id:DUP \\\n")
	shellScriptFile.write("-bedpe bedpe_file:$WORKING_DIR${SAMPLE_ID}.del.100.sorted.bedpe,weight:3,id:DEL \\\n")
	shellScriptFile.write(" -pe  bam_file:$WORKING_DIR$BAM_FILE'.sorted.bam',id:$SAMPLE_ID,histo_file:$WORKING_DIR$BAM_FILE'.histo',read_length:%s,min_non_overlap:100,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:40,mean:$MEAN,stdev:$STDEV \\\n"%args.readlength)
	shellScriptFile.write(" -sr id:$SAMPLE_ID,bam_file:$WORKING_DIR$BAM_FILE'_splitters.sorted.bam',back_distance:10,weight:1,min_mapping_threshold:40 > $WORKING_DIR$SAMPLE_ID'_Lumpy.vcf'\n")

# Lumpy Filter/Annotate Deprecated
def LumpyFilterAnnotate(shellScriptFile):
	shellScriptFile.write("\n#Genotype,Filter,Annotate\n##Run SVtyper\n")
	shellScriptFile.write("/opt/tools/svtyper-0.1.0/svtyper -B $WORKING_DIR$BAM_FILE'.sorted.bam' -i $WORKING_DIR$SAMPLE_ID'_Lumpy.vcf' -o $WORKING_DIR$SAMPLE_ID'_Lumpy_genotyped.vcf'\n")
	shellScriptFile.write("\n##Filter genotyped Lumpy (Custom script)\n")
	shellScriptFile.write("python /opt/tools/VariantAnnotation/FilterLumpyGenotyped.py $WORKING_DIR$SAMPLE_ID'_Lumpy_genotyped.vcf'\\\n")
	shellScriptFile.write(" $WORKING_DIR$SAMPLE_ID'_Lumpy_genotyped_filtered.vcf' $WORKING_DIR$SAMPLE_ID'_Lumpy_DelDup.bed' \\\n")
	shellScriptFile.write(" $WORKING_DIR$SAMPLE_ID'_Lumpy_DelDup.vcf' 0 3000000 \n")
	shellScriptFile.write("\n##Convert to ANNOVAR AVinput\n")
	shellScriptFile.write("perl /opt/tools/annovar/convert2annovar.pl -format vcf4 -allsample --comment --includeinfo -withfreq $WORKING_DIR$SAMPLE_ID'_Lumpy_genotyped_filtered.vcf' > $WORKING_DIR$SAMPLE_ID'_Lumpy.avinput'\n")
	shellScriptFile.write("\n##Run ANNOVAR DGV\n")
	shellScriptFile.write("/mnt/causes-vnx1/DATABASES/annovar/new_table_annovar.pl --otherinfo --buildver hg19 --protocol refgene,dgvMerged_commafix --operation g,r --argument '','-minqueryfrac 0.8 --colsWanted 2&3&4&10&17&18&19&21' \\\n")
	shellScriptFile.write("$WORKING_DIR$SAMPLE_ID'_Lumpy.avinput' /mnt/causes-vnx1/DATABASES/annovar/humandb -out $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar' \n")
	shellScriptFile.write("python /opt/tools/VariantAnnotation/PruneAnnotatedCNVs.py  $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar.hg19_multianno.txt' $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar.hg19_multianno.filtered.txt' \n")
	shellScriptFile.write("python /opt/tools/VariantAnnotation/AddSummaryToAnnovar_CNV.py $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar.hg19_multianno.filtered.txt' $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar.hg19_multianno.filtered.omimAndSummary.txt' \n")

# Not using this anymore
def SVANNOTATE(shellScriptFile):
	shellScriptFile.write("\n# Get variants in proband \n")
	shellScriptFile.write("# Convert to bed using annovar \n")
	shellScriptFile.write("perl /mnt/causes-vnx1/DATABASES/annovar/convert2annovar.pl $WORKDIR/${SAMPLE}.DEL.preformat \\\n")
	shellScriptFile.write("\t -format vcf4 \\\n")
	shellScriptFile.write("\t -outfile $WORKING_DIR/${SAMPLE}.DEL \\\n")
	shellScriptFile.write("\t -includeinfo \n")


# This is the main part of the program. Here I'll parse out what's necessary to run, and add it to the script sequentially
def Main():

	args=GetArgs()

	# Define some variables here based on the arguments you just read in
	sampleID = args.sampleID
	workingDir = args.workingDir
	R1fastq = args.inputR1
	R2fastq = args.inputR2
	mtoolboxConfigFile = args.mtoolbox

	print "The working directory is: %s"%workingDir
	print "The sample ID is: %s"%sampleID
	print "The fastq files you're working with are: \n%s\n%s\n"%(R1fastq,R2fastq)
	
	numProcessors = args.processors
	Memory = args.memory
	runTime = "200:00:00"
	
	
	#######################
	### Write the Script###
	#######################
	
	#This is where we'll begin writing the script, starting with the header information
	shellScriptFile = open('%s%s_FullPipeline_%s.sh'%(workingDir,sampleID,args.version),'w')
	shellScriptFile.write('#!/bin/bash\n')


	##########
	# Header #
	##########
	
	# These will be different for SGE (old cluster) and PBS (new cluster)
	if args.scheduler=='SGE':
		#Job name
		shellScriptFile.write('#$ -N %s_%s_PrimaryPipeline\n'%(sampleID,args.version))
		#Export environment variables
		shellScriptFile.write('#$ -V\n')
		#set location for log files
		shellScriptFile.write('#$ -o %s%s.o\n'%(workingDir,sampleID))
		shellScriptFile.write('#$ -e %s%s.e\n'%(workingDir,sampleID))
		#email on job abort
		shellScriptFile.write('#$ -m bea\n#$ -M prichmond@cmmt.ubc.ca\n')
		#SGE resource requirements
		shellScriptFile.write('#$ -l h_rt=%s,h_vmem=%s\n'%(runTime,Memory))
		#Specific node being submitted to, only if it's specified
		if args.nodeName is not None:
			shellScriptFile.write("#$ -l h=%s\n"%args.nodeName)
		#Parallel processing
		shellScriptFile.write('#$ -pe smp %d\n'%numProcessors)
		shellScriptFile.write('\n\nexport PARALLEL=$NSLOTS\nexport OMP_NUM_THREADS=$NSLOTS\n')
	
	elif args.scheduler == 'PBS':
		#Job name
		shellScriptFile.write('#PBS -N %s_%s_PrimaryPipeline\n'%(sampleID,args.version))
		#Export environment variables
		shellScriptFile.write('#PBS -V\n')
		#set location for log files
		shellScriptFile.write('#PBS -o %s%s%s.o\n'%(workingDir,sampleID,args.version))
		shellScriptFile.write('#PBS -e %s%s%s.e\n'%(workingDir,sampleID,args.version))
		#email on job abort
		if args.Email:
			shellScriptFile.write('#PBS -m bea\n#PBS -M %s\n'%args.Email)
		#SGE resource requirements
		shellScriptFile.write("## Set the total memory for the job\n")
		shellScriptFile.write('#PBS -l mem=%s\n'%Memory)
		shellScriptFile.write("## Set the max walltime for the job\n")
		shellScriptFile.write('#PBS -l walltime=%s\n'%runTime)
		#Specific node being submitted to, only if it's specified
		if args.nodeName is not None:
			shellScriptFile.write("## Not sure we want to specify nodes on new cluster yet\n##PBS -l h=%s\n"%args.nodeName)
		#Parallel processing
		shellScriptFile.write("## Set the total number of processors for the job\n")
		shellScriptFile.write('#PBS -l nodes=1:ppn=%s\n'%numProcessors)
		shellScriptFile.write('NSLOTS=$PBS_NUM_PPN\n')
		shellScriptFile.write('umask 0002\n')
		shellScriptFile.write("source /opt/tools/hpcenv.sh\n\n")
	
	elif args.ascheduler == 'SLURM':
		shellScriptFile.write("#!/bin/bash\n")
		if args.Email:
			shellScriptFile.write("#SBATCH --mail-user=%s\n"%args.Email)
			shellScriptFile.write("#SBATCH --mail-type=ALL\n")
		shellScriptFile.write("#SBATCH --mem=%s\n"%Memory)
		shellScriptFile.write("#SBATCH --cpus-per-task=%s\n"%numProcessors)
		shellScriptFile.write("#SBATCH --time=%s\n"%runTime)
		shellScriptFile.write("#SBATCH --nodes=1\n")
		shellScriptFile.write("#SBATCH --output=%s-%j.out\n")
		shellScriptFile.write("#SBATCH --error=%s-%j.out\n")
	
	#############
	# Variables #
	#############

	#Set the variables for working directory, sampleID, and fastq full filepath
	if args.version=='new':
		shellScriptFile.write("\nSAMPLE_ID=\'%s\'\n"%sampleID)
	
	elif args.version=='old':
		shellScriptFile.write("\nSAMPLE_ID=\'%s_bowtie2\'\n"%sampleID)
		shellScriptFile.write("BOWTIE2_INDEX=\'/mnt/causes-vnx1/GENOMES/hg19/hg19\'\n")
	
	else:
		print "You have specified a wrong value for version.  Must be either old or new"
		print "Unless you're running with 'neither'"
		shellScriptFile.write("\nSAMPLE_ID=\'%s\'\n"%sampleID)
	
	if args.GENOME=='hg19':
	        shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-vnx1/GENOMES/hg19/hg19_bwa.fa\'\n")
	        shellScriptFile.write("BWA_INDEX=\'/mnt/causes-vnx1/GENOMES/hg19/hg19_bwa\'\n")
	        shellScriptFile.write("GENOMEFILE=\'/mnt/causes-vnx1/GENOMES/hg19/hg19_bwa.genome\'\n")
	        shellScriptFile.write("CHROM=\'/mnt/causes-vnx1/GENOMES/hg19/FASTA/\'\n")
	
	elif args.GENOME=='GSC':
	        shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa\'\n")
		shellScriptFile.write("BWA_INDEX=\'/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa\'\n")
		shellScriptFile.write("CHROM=\'/mnt/causes-vnx1/GENOMES/GSC/SplitByChrom/\'\n")
		shellScriptFile.write("GENOMEFILE=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.genome\n\n")
	else:
	        print "You did not choose a viable genome version, choose either GSC or hg19"
	        sys.exit()
	
	shellScriptFile.write("WORKING_DIR=\'%s\'\n"%workingDir)
	shellScriptFile.write("FASTQR1=\'%s\'\n"%R1fastq)
	shellScriptFile.write("FASTQR2=\'%s\'\n"%R2fastq)


	#####################
	# Pipeline Commands #
	#####################
	

	if args.version == 'new':
		shellScriptFile.write("\n echo \"Primary Analysis Started\"\n")
		shellScriptFile.write("date\n")
		FastQC(shellScriptFile)
	        BWA(shellScriptFile)
		Sam2Bam(shellScriptFile)
	        DupRemove(shellScriptFile)
	        Realign(shellScriptFile)
	        CleanUp(shellScriptFile)
	        GATK_HaplotypeCaller(shellScriptFile)
	        GATK_HaplotypeCallerGVCF(shellScriptFile)
	        ValidateSAM(shellScriptFile)
		if args.mtoolbox != "":
			shellScriptFile.write("\necho \"Mitochondrial Variant Analysis Started\"\n")
			shellScriptFile.write("date\n")
			MToolBox(shellScriptFile)
	
		#### 	
		if args.Type == 'Genome'
		MosDepth_WGS(shellScriptFile)
			if args.cnv:
				shellScriptFile.write("\necho \"CNV Analysis Started\"\n")
			        shellScriptFile.write("date\n")
				CNVNATOR(shellScriptFile)
				ERDS(shellScriptFile)

			# Add Smoove?	
			if args.sv:
				if not args.cnv:
					CNVNATOR(shellScriptFile)
				Lumpy(shellScriptFile)
				Pindel(shellScriptFile)
				MetaSV(shellScriptFile)
				shellScriptFile.write("\necho \"CNV Analysis Finished\"\n")
			        shellScriptFile.write("date\n")

		elif args.Type == 'Exome':
			if args.metrics_exome:
				shellScriptFile.write("\necho \"Exome Metrics Calculations Started\"\n")
				shellScriptFile.write("EXOME_CAPTURE_BED=/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/Agilent_SureSelect_Human_All_Exon_V4/S03723314_Covered_chrnameswithoutchr.bed \n")
				shellScriptFile.write("EXOME_CAPTURE_INTERVAL=/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/Agilent_SureSelect_Human_All_Exon_V4/S03723314_Covered_chrnameswithoutchr.GRCh37-lite.interval_list \n")
				shellScriptFile.write("METRICS_WORKING_DIR=$WORKING_DIR/METRICS/ \n")
				shellScriptFile.write("mkdir -p $METRICS_WORKING_DIR \n")
				shellScriptFile.write(" \n")
	
				MosDepth_WES(shellScriptFile)
				Picard_HSMETRICS(shellScriptFile)
	
	
		shellScriptFile.write("\n echo \"Primary Analysis Finished\"\n")
		shellScriptFile.write("date\n")
	
	# Not changing anything for exome vs. genome here	
	elif args.version == 'old':
		FastQC(shellScriptFile)
		Bowtie2(shellScriptFile)
		Sam2Bam(shellScriptFile)
		DupRemove(shellScriptFile)
	        Realign(shellScriptFile)
	        VarCall(shellScriptFile)
	        SummaryStats(shellScriptFile)
	        IndelRemoval(shellScriptFile)
		ReadMetrics(shellScriptFile)
		SNVMetrics(shellScriptFile)
		CleanUp(shellScriptFile)
		ValidateSAM(shellScriptFile)
		GATK_Coverage(shellScriptFile)
	
	elif args.version == 'neither':
		print "Not putting down the primary analysis"
	
	


if __name__=="__main__":
	Main()
