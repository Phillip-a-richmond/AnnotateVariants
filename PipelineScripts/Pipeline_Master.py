import sys, os, argparse

#####################################
# Author:   Phillip Richmond        #
# Contact:  prichmond@cmmt.ubc.ca   #
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
# March 25th 2019
	# Removed SV Calling from this pipeline. Scripts getting too dirty and long
	# SV Calling is now handled within Bam2SV.py	

# September 9th 2020
        # Overhaul of the code for the new GPCC cluster.
        # This means: 
            # Adding SLURM
            # Moving filepaths
            # changing assumed DATABASES directory
            # Dropping hg19, adding GRCh38-lite
            


##################
### Initialize ###
##################


# Get the args from the user
def GetArgs():
	if len(sys.argv) < 2:
		print "Re-run with the -h option"
		print "Typical Running Command (genome singleton with mito analysis):"
		print "python Pipeline_Master.py -p 16 -m 60G -R 150 -T Genome -S PBS --mtoolbox /mnt/causes-vnx1/PIPELINES/AnnotateVariants/MToolBox_config_files/MToolBox_rCRS_config_with_markdup_and_indelrealign_RvdL.sh -v new  -s TIDEX555 -d /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T149/  -1 /mnt/causes-vnx2/TIDE/RAW/GENOME_TIDEX/T149/TIDEX555_R1.fastq.gz  -2 /mnt/causes-vnx2/TIDE/RAW/GENOME_TIDEX/T149/TIDEX555_R2.fastq.gz -G GSC"

		print "Typical Running Command (exome trio)"
		print "python Pipeline_Master.py -p 16 -m 60G -R 150 -T Exome -G GSC -S PBS -s TIDEX698 -v new -d /mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/T186/ -1 /mnt/causes-vnx2/TIDE/RAW/EXOME_TIDEX/T186/TIDEX698_1.fastq.gz -2  /mnt/causes-vnx2/TIDE/RAW/EXOME_TIDEX/T186/TIDEX698_2.fastq.gz"
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
	parser.add_argument("-R","--readlength",help="The length of the reads in integer form",default=100)
	parser.add_argument("-T","--Type",help="Exome || Genome",required=True)
	parser.add_argument("-S","--scheduler",help="Which scheduler you want to submit to.  This will determine the format of the shell script. Options: PBS || SGE",type=str)
	parser.add_argument("-G","--GENOME",help="Which Genome version do you want to use? Options are GSC || hg19",required=True)
	parser.add_argument("-V","--VCF",help="Will run GATK Variant calling to produce a VCF",action='store_true')
	parser.add_argument("--mtoolbox",help="Provide an MToolBox config file here to perform a mitochondrial variant analysis, e.g. /path/to/AnnotateVariants/MToolBox_config_files/MToolBox_rCRS_config_with_markdup_and_indelrealign_RvdL.sh",default="/mnt/causes-vnx1/PIPELINES/AnnotateVariants/MToolBox_config_files/MToolBox_rCRS_config_with_markdup_and_indelrealign_RvdL.sh")
	parser.add_argument("--metrics_exome",help="Calculate exome coverage and other metrics using Mosdepth and Picard CalculateHsMetrics, assuming the Agilent_SureSelect_Human_All_Exon_V4 capture kit was used",action='store_true',default=False)
	parser.add_argument("--metrics_genome",help="Calculate genome coverage and other metrics using Mosdepth and Picard CalculateHsMetrics",action='store_true',default=False)
	parser.add_argument("-E","--Email",help="Email address",type=str)
	args = parser.parse_args()
	return args



#########################################################################################


# Now define these functions, which when called, they will write out the command to the shell
# script.  This is an easy way to control what goes into the script, so if I have to run
# only fastQC on all the files, then I can just specify that, and call that command

###### Either new pipeline, or shared pipeline commands (both old and new) #######

# FastQC to generate quality report. No interpretation of the quality report though
def FastQC(shellScriptFile):
	shellScriptFile.write("\n#FastQC\n")
	shellScriptFile.write("$FASTQC --extract $FASTQR1 $FASTQR2 -o $WORKING_DIR\n")

#Map with BWA
def BWA(shellScriptFile):
	shellScriptFile.write("\n#Map with BWA\n")
	shellScriptFile.write("$BWA mem $BWA_INDEX -t $NSLOTS -R \"@RG\\tID:$SAMPLE_ID\\tSM:$SAMPLE_ID\\tPL:illumina\" -M $FASTQR1 $FASTQR2 > $WORKING_DIR$SAMPLE_ID'.sam'\n\n")

#Convert to binary, sort, and index
def Sam2Bam(shellScriptFile):
	shellScriptFile.write("#Convert to binary, sort, and index\n")
	shellScriptFile.write("$SAMTOOLS view -@ $NSLOTS -u -bS $WORKING_DIR$SAMPLE_ID\'.sam\' | $SAMTOOLS sort -@ $NSLOTS -m 3G - $WORKING_DIR$SAMPLE_ID\'.sorted\'\n")
	shellScriptFile.write("$SAMTOOLS index $WORKING_DIR$SAMPLE_ID\'.sorted.bam\'\n")

# Duplicate marking
def DupRemove(shellScriptFile):
	shellScriptFile.write("\n#Remove Duplicates\n")
	shellScriptFile.write("TMPDIR=$WORKING_DIR'picardtmp/'\n")
	shellScriptFile.write("mkdir $TMPDIR\n")
	shellScriptFile.write("$JAVA -jar $PICARD MarkDuplicates I=$WORKING_DIR$SAMPLE_ID\'.sorted.bam\' O=$WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\' REMOVE_DUPLICATES=false TMP_DIR=$TMPDIR M=$WORKING_DIR$SAMPLE_ID\'_DuplicateResults.txt\'\n")
	shellScriptFile.write("$SAMTOOLS index $WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\'\n")

# GATK Realignment
def Realign(shellScriptFile):
	shellScriptFile.write("\n#Realign\n")
	shellScriptFile.write("$JAVA -jar $GATKJAR -T RealignerTargetCreator -nt $NSLOTS -R $GENOME_FASTA -minReads 5 -I $WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_indelsites.intervals\' \n")
	shellScriptFile.write("$JAVA -jar $GATKJAR -T IndelRealigner -model USE_READS -R $GENOME_FASTA -targetIntervals $WORKING_DIR$SAMPLE_ID\'_indelsites.intervals\' -I $WORKING_DIR$SAMPLE_ID\'_dupremoved.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\'\n")
	shellScriptFile.write("$SAMTOOLS index $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\'\n")

# Added DJA 2016/06/15
def ReadMetrics(shellScriptFile):
	shellScriptFile.write("\n#Read Metrics\n")
	shellScriptFile.write("$JAVA -jar $PICARD CollectMultipleMetrics INPUT=$WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\' OUTPUT=$WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted\' REFERENCE_SEQUENCE=$GENOME_FASTA\n")

# Validate SAM file
def ValidateSAM(shellScriptFile):
	shellScriptFile.write("\n#Picard Validate BAM\n")
	shellScriptFile.write("$JAVA -jar $PICARD ValidateSamFile M=VERBOSE MO=10000000 I=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' O=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam_PicardSAMValidate.txt'\n")

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
	shellScriptFile.write("$JAVA -jar $GATKJAR  -T HaplotypeCaller \\\n")
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
	shellScriptFile.write("/opt/tools/htslib-1.5/bin/tabix $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.vcf.gz \n")
	shellScriptFile.write("\n#Filter\n")
	shellScriptFile.write("/opt/tools/bcftools-1.5/bin/bcftools filter \\\n")
	shellScriptFile.write("--include 'FORMAT/AD[1]>=3 && INFO/MQ>=20 && %QUAL>=20' \\\n")
	shellScriptFile.write("$WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.vcf.gz \\\n")
	shellScriptFile.write("--output-type z \\\n")
	shellScriptFile.write("--output $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.norm.filter.vcf.gz \n")

#SNP calling, creating a GVCF which is used in multi-sample variant calling (Bam2Gemini.sh part of the pipeline)
def GATK_HaplotypeCallerGVCF(shellScriptFile):
	shellScriptFile.write("\n#SNP Calling HaplotypeCaller GVCFmode\n")
	shellScriptFile.write("$JAVA -jar $GATKJAR \\\n")
	shellScriptFile.write(" -T HaplotypeCaller --emitRefConfidence GVCF \\\n")
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
	shellScriptFile.write("$JAVA -jar $GATKJAR -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $GENOME_FASTA -I $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_dupremoved_realigned_Coverage'\n")

# MosDepth depth of coverage assessment
def MosDepth_WGS(shellScriptFile,args):
	shellScriptFile.write("\n#MosDepth \n\n")
	shellScriptFile.write("$MOSDEPTH -n -q 0:1:10:20:40:60:100: -t $NSLOTS \\\n")
	shellScriptFile.write("$WORKING_DIR$SAMPLE_ID $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	shellScriptFile.write("\n#Plot output \n\n")
	shellScriptFile.write("python /opt/tools/mosdepth-0.2.2/plot-dist.py \\\n")
	shellScriptFile.write("$WORKING_DIR${SAMPLE_ID}.mosdepth.global.dist.txt \\\n")
	shellScriptFile.write("-o $WORKING_DIR${SAMPLE_ID}.mostdepthCoverage.html \\\n")
	shellScriptFile.write("> $WORKING_DIR${SAMPLE_ID}.mostdepthCoverage.summary.txt \n\n")

def MosDepth_WES(shellScriptFile,args):
	shellScriptFile.write("\n# Run MosDepth \n")
	shellScriptFile.write("SAMPLE=\'%s\'\n"%args.sampleID)
	shellScriptFile.write("EXOME_CAPTURE_BED=/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/Agilent_SureSelect_Human_All_Exon_V4/S03723314_Covered_chrnameswithoutchr.bed \n")
	shellScriptFile.write("MOSDEPTH_PATH=/opt/tools/mosdepth-0.2.2/ \n")
	shellScriptFile.write(" \n")
	shellScriptFile.write("$MOSDEPTH_PATH/mosdepth -t 4 -n -b $EXOME_CAPTURE_BED $METRICS_WORKING_DIR/$SAMPLE $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	shellScriptFile.write("python $MOSDEPTH_PATH/plot-dist.py -o $METRICS_WORKING_DIR/mosdepth_coverage.html $METRICS_WORKING_DIR/*region.dist.txt \n")
	shellScriptFile.write("COV_OUT=$METRICS_WORKING_DIR/mean_coverage_mosdepth.txt \n")
	shellScriptFile.write("python $MOSDEPTH_PATH/mean_coverage_mosdepth.py $METRICS_WORKING_DIR/*regions.bed.gz > $COV_OUT \n")
	shellScriptFile.write("head -1000 $COV_OUT \n")
	shellScriptFile.write(" \n")

def Picard_HSMETRICS_genome(shellScriptFile):
	shellScriptFile.write("\n# Run Picard CalculateHsMetrics \n")
	shellScriptFile.write("$JAVA -jar $PICARD CalculateHsMetrics \\\n")
	shellScriptFile.write("\t I=$WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
	shellScriptFile.write("\t O=$METRICS_WORKING_DIR${SAMPLE_ID}_picard_hsmetrics.txt \\\n")
	shellScriptFile.write("\t R=$GENOME_FASTA \\\n")
	shellScriptFile.write(" \n")

def Picard_HSMETRICS_exome(shellScriptFile):
        shellScriptFile.write("\n# Run Picard CalculateHsMetrics \n")
        shellScriptFile.write("$JAVA -jar $PICARD CalculateHsMetrics \\\n")
        shellScriptFile.write("\t I=$WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
        shellScriptFile.write("\t O=$METRICS_WORKING_DIR${SAMPLE_ID}_picard_hsmetrics.txt \\\n")
        shellScriptFile.write("\t R=$GENOME_FASTA \\\n")
        shellScriptFile.write("\t BAIT_INTERVALS=$EXOME_CAPTURE_INTERVAL \\\n")
        shellScriptFile.write("\t TARGET_INTERVALS=$EXOME_CAPTURE_INTERVAL \n")
        shellScriptFile.write(" \n")

# MToolBox Mitochondrial analysis
def MToolBox(shellScriptFile,sampleID,mtoolboxConfigFile):
	shellScriptFile.write("\necho \"Mitochondrial Variant Analysis Started\"\n")
	shellScriptFile.write("date\n")
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


###### Old Pipeline, Legacy Commands ######

# Part of the old pipeline, maintaining for legacy purposes. Not updating this
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
	shellScriptFile.write("WORKING_DIR=\'%s\'\n"%workingDir)
	shellScriptFile.write("FASTQR1=\'%s\'\n"%R1fastq)
	shellScriptFile.write("FASTQR2=\'%s\'\n"%R2fastq)

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
	
	if args.metrics_exome:
		shellScriptFile.write("EXOME_CAPTURE_BED=/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/Agilent_SureSelect_Human_All_Exon_V4/S03723314_Covered_chrnameswithoutchr.bed \n")
		shellScriptFile.write("EXOME_CAPTURE_INTERVAL=/mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/Agilent_SureSelect_Human_All_Exon_V4/S03723314_Covered_chrnameswithoutchr.GRCh37-lite.interval_list \n")
		shellScriptFile.write("METRICS_WORKING_DIR=$WORKING_DIR/METRICS/ \n")
		shellScriptFile.write("mkdir -p $METRICS_WORKING_DIR \n")
		shellScriptFile.write(" \n")
	


        # Make some variables for tools which are used within this pipeline
        shellScriptFile.write("# Define Tool paths. If they are in your path, simply change these full filepaths to only be the final command\n")
        shellScriptFile.write("# For example: Change BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools to be BCFTOOLS=bcftools if it's in your path \n\n")
        shellScriptFile.write("BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools\n")
        shellScriptFile.write("GATKJAR=/opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar\n")
        shellScriptFile.write("JAVA=/opt/tools/jdk1.7.0_79/bin/java\n")
        shellScriptFile.write("BGZIP=/opt/tools/tabix/bgzip\n")
        shellScriptFile.write("TABIX=/opt/tools/tabix/tabix\n")
	shellScriptFile.write("BWA=/opt/tools/bwa-0.7.12/bwa\n")
	shellScriptFile.write("SAMTOOLS=/opt/tools/samtools-1.2/samtools\n")
	shellScriptFile.write("PICARD=/opt/tools/picard-tools-1.139/picard.jar\n")
	shellScriptFile.write("MOSDEPTH=/opt/tools/mosdepth-0.2.2/mosdepth\n")
	shellScriptFile.write("FASTQC=/opt/tools/FastQC/fastqc\n")

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
		if args.VCF:
	        	GATK_HaplotypeCaller(shellScriptFile)
	        GATK_HaplotypeCallerGVCF(shellScriptFile)
	        ValidateSAM(shellScriptFile)
		if args.mtoolbox != "":
			MToolBox(shellScriptFile,sampleID,mtoolboxConfigFile)
		if (args.Type == "Genome"):
			if args.metrics_genome:
				MosDepth_WGS(shellScriptFile,args)
				Picard_HSMETRICS_genome(shellScriptFile)
		elif args.Type == 'Exome':
			if args.metrics_exome:
				shellScriptFile.write("\necho \"Exome Metrics Calculations Started, assuming the Agilent_SureSelect_Human_All_Exon_V4 capture kit\"\n")
				MosDepth_WES(shellScriptFile,args)
				Picard_HSMETRICS_exome(shellScriptFile)
	
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
