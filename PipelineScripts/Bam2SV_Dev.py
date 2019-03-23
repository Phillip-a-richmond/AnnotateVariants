import sys, os, argparse

#####################################
# Author:   Phillip Richmond        #
# Contact:  prichmond@cmmt.ubc.ca   #
# Open source GNU licensing         #
#####################################

########
#Bam2SV#
########

# The purpose of this program is to take in a BAM or set of BAM files, and call variants that aren't SNVs and indels
# It will support:
	# Copy Number Variants - CNVnator / ERDS
	# Structural Variants - SMOOVE
	# Mobile Element Insertions - MELT
	# Short tandem repeats - STRetch
	# Mitochondiral Analysis - MToolbox 

# Initially, I'm just going to implement the singleton approach. Joint genotyping, or merging genotypes over these variants is a little more complex

##################
### Initialize ###
##################
def GetArgs():
	if len(sys.argv) < 2:
		print "Typical Running Command: Singleton BAM file: CNV, SV, MEI, STR"
	        print "python Bam2SV.py  -d /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T274_new/ -p 16 -m 40G \\"
	        print "-B T014-1_BWAmem_dupremoved_realigned.sorted.bam\n"
		sys.exit()
	
	# Read in your arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-d","--args.workingDir",help="Working directory on your filesystem",required=True)
	parser.add_argument("-F","--Family",help="The Family Identifier",required=True)
	parser.add_argument("-p","--processors",help="Choose the number of processors for this job",type=int,required=True)
	parser.add_argument("-T","--Type",help="Genome || Exome",required=True)
	parser.add_argument("-m","--memory",help="Choose the memory needed for this job",required=True)
	parser.add_argument("-B","--BAMLIST",help="Comma separated list of BAMS to be processed.  These BAMs should correspond to the identifiers inside the Ped File")
	parser.add_argument("-E","--Email",help="Email address",type=str)
	parser.add_argument("-G","--GENOME",help="Which Genome version do you want to use? Options are GSC || hg19",required=True)
	parser.add_argument("-S","--Singleton",help="Use this option if you are running a singleton",action='store_true',default=False)
	args = parser.parse_args()
	return args	


#########################################################################################

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


# Mobile element insertion calling
def MEI(shellScriptFile):
        shellScriptFile.write("\n# Mobile Element Insertions\n")
        shellScriptFile.write("## Define Variables\n")
        shellScriptFile.write("ANALYSIS_DIR=${WORKING_DIR}MEI\n")
        shellScriptFile.write("MELT_DIR=/opt/tools/MELVTv2.1.5/\n")
        shellScriptFile.write("MEI_LIST=${MELT_DIR}/me_refs/1KGP_Hg19/mei_list.txt\n")
        shellScriptFile.write("GENE_ANNO=/opt/tools/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed\n\n")
        shellScriptFile.write("# MELT Singleton \n\n")
        shellScriptFile.write("java -jar ${MELT_DIR}MELT.jar Single \\\n")
        shellScriptFile.write("-a -b hs37d5/NC007605 -c 8 -h $GENOME_FASTA \\\n")
        shellScriptFile.write("-bamfile $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
        shellScriptFile.write("-n ${MELT_DIR}add_bed_files/1KGP_Hg19/hg19.genes.bed \\\n")
        shellScriptFile.write(" -t $MEI_LIST -w $WORKING_DIR \n")
# This is for multiple samples, leaving here for now
        shellScriptFile.write("# MELT Assuming multiple samples downstream \n")
        shellScriptFile.write("## Step 1 - Preprocess \n")
        shellScriptFile.write("java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \\\n")
        shellScriptFile.write("-bamfile $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
        shellScriptFile.write("-h $GENOME_FASTA \n\n")
        shellScriptFile.write("## Step 2 - Individual Analysis\n")
        shellScriptFile.write("java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \\\n")
        shellScriptFile.write("-w $ANALYSIS_DIR \\\n")
        shellScriptFile.write("-t $MEI_LIST -c 30 -h $GENOME_FASTA \\\n")
        shellScriptFile.write("-bamfile $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")


def RunSTRetch(shellScriptFile):
	shellScriptFile.write("# Running STRetch")
	shellScriptFile.write("/mnt/causes-vnx1/PIPELINES/STRetch/tools/bin/bpipe run \\\n")
	shellScriptFile.write("-p input_regions=/mnt/causes-vnx1/PIPELINES/STRetch/reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed \\\n")
	shellScriptFile.write("/mnt/causes-vnx1/PIPELINES/STRetch/pipelines/STRetch_wgs_bam_pipeline.groovy \\\n")
	shellScriptFile.write("$WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	
# This is basic, just running the VCFAnno command on your input merged VCF
def RunVCFAnno(shellScriptFile,DatabaseDir):
	shellScriptFile.write("\n# Step 5: VCFAnno - Turn your VCF file into an annotated VCF file\n")
	shellScriptFile.write('/opt/tools/vcfanno/vcfanno -lua /mnt/causes-vnx1/PIPELINES/AnnotateVariants/VCFAnno/custom.lua \\\n')
	shellScriptFile.write('-p $NSLOTS -base-path %s\\\n'%DatabaseDir)
	shellScriptFile.write('/mnt/causes-vnx1/PIPELINES/AnnotateVariants/VCFAnno/VCFAnno_Config_20190131_GAC.toml \\\n')
	shellScriptFile.write('$NORMFILTERVCF > $ANNOVCF \n\n')

# NOTE: If you want to add certain things as --a-ok make sure you add them here, otherwise they may error on the creation of the mysqlDB
def VCF2DB(shellScriptFile):
	shellScriptFile.write("\n# Step 6: VCF2DB - Turn your annotated VCF file into a GEMINI DB\n\n")
	shellScriptFile.write('python /opt/tools/vcf2db/vcf2db.py \\\n')
	shellScriptFile.write('--expand gt_quals --expand gt_depths --expand gt_alt_depths --expand gt_ref_depths --expand gt_types \\\n')
	shellScriptFile.write(' --a-ok InHouseDB_AC  --a-ok in_segdup --a-ok AF --a-ok AC --a-ok AN --a-ok MLEAC --a-ok MLEAF --a-ok gnomad_genome_hom_global --a-ok gnomad_genome_hom_afr --a-ok gnomad_genome_hom_amr --a-ok gnomad_genome_hom_asj --a-ok gnomad_genome_hom_eas --a-ok gnomad_genome_hom_fin --a-ok gnomad_genome_hom_nfe --a-ok gnomad_genome_hom_oth --a-ok gnomad_exome_hom_global --a-ok gnomad_exome_hom_afr --a-ok gnomad_exome_hom_amr --a-ok gnomad_exome_hom_asj --a-ok gnomad_exome_hom_eas --a-ok gnomad_exome_hom_fin --a-ok gnomad_exome_hom_nfe --a-ok gnomad_exome_hom_oth --a-ok cpg_island --a-ok common_pathogenic --a-ok cse-hiseq --a-ok DS --a-ok ConfidentRegion \\\n')
	shellScriptFile.write('$ANNOVCF $PED_FILE $GEMINIDB \n')

def Singleton_RenameVCF2MergedVCF(shellScriptFile):
	shellScriptFile.write("\n# Step 1-2: rename your VCF to merged VCF\n")
	shellScriptFile.write("cp $WORKING_DIR$SAMPLE1_VCF $WORKING_DIR${FAMILY_ID}.merged.hc.vcf\n")

def Main():

	args=GetArgs()
	print args
#######################
### Write the Script###
#######################

	#This is where we'll begin writing the script, starting with the header information
	shellScriptFile = open('%s%s_Bam2Gemini.sh'%(args.workingDir,args.Family),'w')
	shellScriptFile.write('#!/bin/bash\n')
	
	# These are parameters for the scheduler
	numProcessors = args.processors
	Memory = args.memory
	runTime = "240:00:00"

	#Job name
	shellScriptFile.write('#PBS -N %s_Bam2Gemini\n'%args.Family)
	#Export environment variables
	shellScriptFile.write('#PBS -V\n')
	#set location for log files
	shellScriptFile.write('#PBS -o %s%s.o\n'%(args.workingDir,args.Family))
	shellScriptFile.write('#PBS -e %s%s.e\n'%(args.workingDir,args.Family))
	#email on job abort
	if args.Email:
		shellScriptFile.write('#PBS -m bea\n#PBS -M %s\n'%args.Email)
	shellScriptFile.write("## Set the total memory for the job\n")
	shellScriptFile.write('#PBS -l mem=%s\n'%Memory)
	shellScriptFile.write("## Set the max walltime for the job\n")
	shellScriptFile.write('#PBS -l walltime=%s\n'%runTime)
	shellScriptFile.write("## Set the total number of processors for the job\n")
	shellScriptFile.write('#PBS -l nodes=1:ppn=%s\n'%numProcessors)
	shellScriptFile.write('NSLOTS=$PBS_NUM_PPN\n')
	shellScriptFile.write('umask 0002\n')
	shellScriptFile.write('source /opt/tools/hpcenv.sh\n\n')
	
	
	
		
	#Set the variables for working directory, sampleID, and fastq full filepath
	shellScriptFile.write("FAMILY_ID=\'%s\'\n"%args.Family)
	shellScriptFile.write("WORKING_DIR=\'%s\'\n"%args.workingDir)
	if args.GENOME=='hg19':
		shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa\'\n")
	elif args.GENOME=='GSC':
		shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa\'\n")
	else:
		print "You did not choose a viable genome version, choose either GSC or hg19"
		sys.exit(shellScriptFile)
	shellScriptFile.write("PED_FILE=$WORKING_DIR/%s\n"%args.PED)
	shellScriptFile.write("TMPDIR=${WORKING_DIR}tmpdir/\n")
	shellScriptFile.write("mkdir $TMPDIR\n")
	
	#The BAM files from our list, if running in the BAM mode
	if args.BAMLIST:
		for i in range(1,len(BAMS)+1,1):
			shellScriptFile.write("SAMPLE%d_BAM=%s\n"%(i,BAMS[i-1]))
	
	# The VCF files from our list, if running in the VCF mode
	if args.VCFLIST:
		for i in range(1,len(VCFS)+1,1):
			shellScriptFile.write("SAMPLE%d_VCF=%s\n"%(i,VCFS[i-1]))
		# Make sure you have the correct working Directory.
		args.workingDir = args.args.workingDir
		print "The working directory is: %s"%args.workingDir
		print "Family you're working with: %s"%args.Family
		
		if args.BAMLIST:
			BAMS = args.BAMLIST.split(',')
			print "These are the BAM files you're working with:"
			for each in BAMS:
				print each
		elif args.VCFLIST:
			VCFS = args.VCFLIST.split(',')
		        print "These are the VCF files you're working with:"
		        for each in VCFS:
		                print each
			if args.vcftype == 'GVCF':
				print "These should be GVCFs"
			elif args.vcftype == 'VCF':
				print "These should be VCFs"
			else:
				print "You're not using the proper vcftype, or you haven't set it."
				sys.exit(shellScriptFile)
	
	# Add variable names to the shell script
	shellScriptFile.write("# Define some variables\n\n")
	shellScriptFile.write("SNPEFFJAR=/opt/tools/snpEff/snpEff.jar\n")
	shellScriptFile.write("GEMINIDB=$WORKING_DIR${FAMILY_ID}.db\n")
	shellScriptFile.write("VCF=$WORKING_DIR${FAMILY_ID}.merged.hc.vcf\n")
	shellScriptFile.write("NORMVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcf.gz\n")
	shellScriptFile.write("NORMFILTERVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.filter.vcf.gz\n")
	shellScriptFile.write('ANNOVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcfanno.vcf.gz \n')
	

#### Call functions to populate the commands in the script

	# Run in either VCF or BAM mode
	if args.VCFLIST:
	# if running in GVCF mode, then you want to joint genotype from multiple GVCFs
		if args.vcftype == 'GVCF':
			if args.Singleton:
				print "You cannot use singleton with GVCF mode"
				sys.exit(shellScriptFile)
			MergeGVCF_withVCFLIST(shellScriptFile)
			MergedVCF2NormVCF(shellScriptFile)
			FilterVCF(shellScriptFile)
			RunVCFAnno(shellScriptFile,DatabaseDir)
			VCF2DB(shellScriptFile)
	# if running in the VCF mode, I assume you have separate VCFs whcih you want to merge into a single merged.vcf
		elif args.vcftype == 'VCF':
			if args.Singleton(shellScriptFile):
				SingletonVCF2MergedVCFRename(shellScriptFile)
			else:
				MergeVCF_withVCFLIST(shellScriptFile)
			MergedVCF2NormVCF(shellScriptFile)
			FilterVCF(shellScriptFile)
			RunVCFAnno(shellScriptFile,DatabaseDir)
			VCF2DB(shellScriptFile)
	
	elif args.BAMLIST:
		if args.vcftype == 'GVCF':
			if args.Singleton:
                                print "You cannot use singleton with GVCF mode"
                                sys.exit(shellScriptFile)
			Bam2GVCF(shellScriptFile)
			MergeGVCF_withBAMLIST(shellScriptFile)
			MergedVCF2NormVCF(shellScriptFile)
			FilterVCF(shellScriptFile)
			RunVCFAnno(shellScriptFile,DatabaseDir)
			VCF2DB(shellScriptFile)
	
		elif args.vcftype == 'VCF':
			Bam2VCF(shellScriptFile)
			if args.Singleton:
				SingletonVCF2MergedVCFRename(shellScriptFile)
			else:
				MergeVCF_withVCFLIST(shellScriptFile)
	                MergedVCF2NormVCF(shellScriptFile)
	                FilterVCF(shellScriptFile)
			RunVCFAnno(shellScriptFile,DatabaseDir)
			VCF2DB(shellScriptFile)
	
