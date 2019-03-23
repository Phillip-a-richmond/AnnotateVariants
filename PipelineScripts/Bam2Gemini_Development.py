import sys, os, argparse

#####################################
# Author:   Phillip Richmond        #
# Contact:  prichmond@cmmt.ubc.ca   #
# Open source GNU licensing         #
#####################################

############
#Bam2Gemini#
############

# Update January 31st 2019 - Major 2019 overhaul
	# Removed deprecated code
	# Updated VCFAnno config files
	# Updated BCFTools filters
<<<<<<< HEAD
# Update March 12th 2019
	# added singleton
	# added databaseDir, which will have to be fully organized for this release
		# note, databaseDir option will set your database directory, which VCFAnno refers to as it's base-path when looking for annotation files. 
		# This way, annotation files will lie within a single directory
=======
# Update March 1st 2019
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3

##################
### Initialize ###
##################
def GetArgs():
	if len(sys.argv) < 2:
		print "Re-run with the -h option"
		print "Typical Running Command (From Bams, GVCF mode):"
<<<<<<< HEAD
		print "python Bam2Gemini.py  -d /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T149/ -p 16 -m 40G \\"
		print "-P T149.ped -F T149 -v GVCF -G GSC -T Genome \\"
		print "-B TIDEX555_BWAmem_dupremoved_realigned.sorted.bam\n"
	
	
		print "Typical Running Command (From GVCFs):"
	        print "python Bam2Gemini.py  -d /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T149_new/ -p 16 -m 40G \\"
	        print "-P T149.ped -F T149 -v GVCF \\"
	        print "-V T149-1_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf,T149-2_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf,T149-3_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf\n"
=======
		print "python Bam2Gemini.py  -d /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T014_new/ -p 16 -m 40G \\"
		print "-P T014.ped -F T014 -v GVCF \\"
		print "-B T014-1_BWAmem_dupremoved_realigned.sorted.bam,T014-2_BWAmem_dupremoved_realigned.sorted.bam,T014-3_BWAmem_dupremoved_realigned.sorted.bam\n"
	
	
		print "Typical Running Command (From GVCFs):"
	        print "python Bam2Gemini.py  -d /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T014_new/ -p 16 -m 40G \\"
	        print "-P T014.ped -F T014 -v GVCF \\"
	        print "-V T014-1_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf,T014-2_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf,T014-3_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf\n"
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	
		print "Typical Running Command, duo (From VCFs):"
	        print "python Bam2Gemini.py  -d /mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T274_new/ -p 16 -m 40G \\"
	        print "-P T274.ped -F T274 -v VCF \\"
	        print "-V tidex1006_BWAmem_dupremoved_realigned_HaplotypeCaller.vcf,tidex1007_BWAmem_dupremoved_realigned_HaplotypeCaller.vcf\n"
	
		sys.exit()
	
	# Read in your arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-d","--workingDir",help="Working directory on your filesystem",required=True)
	parser.add_argument("-F","--Family",help="The Family Identifier",required=True)
	parser.add_argument("-p","--processors",help="Choose the number of processors for this job",type=int,required=True)
	parser.add_argument("-T","--Type",help="Genome || Exome",required=True)
	parser.add_argument("-m","--memory",help="Choose the memory needed for this job",required=True)
	parser.add_argument("-P","--PED",help="The PED File for this family",required=True)
	parser.add_argument("-B","--BAMLIST",help="Comma separated list of BAMS to be processed.  These BAMs should correspond to the identifiers inside the Ped File")
	parser.add_argument("-E","--Email",help="Email address",type=str)
	parser.add_argument("-v","--vcftype",help="The type of VCF.  If these are GVCFs that need to be merged, or separately called VCFs. REQUIRED option.",required=True)
	parser.add_argument("-V","--VCFLIST",help="Comma separated list of VCF files.  Set either GVCF or VCF with -v to know merging option.")
	parser.add_argument("-G","--GENOME",help="Which Genome version do you want to use? Options are GSC || hg19",required=True)
<<<<<<< HEAD
	parser.add_argument("-S","--Singleton",help="Use this option if you are running a singleton",action='store_true',default=False)
	parser.add_argument("-D","--DatabaseDir",help="Path for the directory containing databases which VCFAnno will use to annotate your VCF. Default is: /mnt/causes-vnx1/DATABASES/",type=str,default='/mnt/causes-vnx1/DATABASES/')
	parser.add_argument("-A","--AnnotateVariantsDir",help="Path for the github repo AnnotateVariants",required=True,default='/mnt/causes-vnx1/PIPELINES/AnnotateVariants')
=======
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	#parser.add_argument("-C","--Config",help="Config file which points at locations for tool executables",required=True)
	#parser.add_argument("-Q","--QueryScript",help="Query script template",required=True)
	args = parser.parse_args()
	return args	


<<<<<<< HEAD
#########################################################################################

# Generate GVCFs from BAM files
def Bam2GVCF(shellScriptFile,BAMS):
	shellScriptFile.write("\n echo \"Primary Analysis Started\"\n")
	shellScriptFile.write("# Generate gVCFs\n\n")
	# Here I'm going to write out a Bam2gVCF for each of the bam files listed.  I'll just add a HC.g.vcf tag to the bam file names
	for i in range(1,len(BAMS)+1,1):
	        shellScriptFile.write("$JAVA -jar $GATKJAR \\\n")
=======



#########################################################################################

# Generate GVCFs from BAM files
def Bam2GVCF():
	shellScriptFile.write("\n echo \"Primary Analysis Started\"\n")
	shellScriptFile.write("# Step 1: Generate gVCFs\n\n")
	# Here I'm going to write out a Bam2gVCF for each of the bam files listed.  I'll just add a HC.g.vcf tag to the bam file names
	for i in range(1,len(BAMS)+1,1):
	        shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar \\\n")
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	        shellScriptFile.write(" -T HaplotypeCaller --emitRefConfidence GVCF \\\n")
		shellScriptFile.write(" -R $GENOME_FASTA -I $WORKING_DIR$SAMPLE%d_BAM \\\n"%(i))
	        shellScriptFile.write(" --standard_min_confidence_threshold_for_calling 10 \\\n")
	        shellScriptFile.write(" --standard_min_confidence_threshold_for_emitting 10 \\\n")
	        shellScriptFile.write(" --min_mapping_quality_score 0 \\\n")
	        shellScriptFile.write(" --min_base_quality_score 10 \\\n")
	        shellScriptFile.write(" --minReadsPerAlignmentStart 5 \\\n")
	        shellScriptFile.write(" --minPruning 2 \\\n")
	        shellScriptFile.write(" --pcr_indel_model NONE \\\n")
	        shellScriptFile.write(" --dbsnp /opt/tools/GATK-3.5-0-g36282e4/resources/dbsnp_138.b37.excluding_sites_after_129.vcf \\\n")
		shellScriptFile.write(" -o $WORKING_DIR${SAMPLE%d_BAM}.HC.g.vcf \n"%i)
# Generate VCFs from BAM files
<<<<<<< HEAD
def Bam2VCF(shellScriptFile,BAMS):
	shellScriptFile.write("\n echo \"Primary Analysis Started\"\n")
	shellScriptFile.write("# Generate VCFs\n\n")
	# Here I'm going to write out a Bam2gVCF for each of the bam files listed.  I'll just add a HC.g.vcf tag to the bam file names
	for i in range(1,len(BAMS)+1,1):
		shellScriptFile.write("$JAVA -jar $GATKJAR -T HaplotypeCaller  \\\n")
=======
def Bam2VCF():
	shellScriptFile.write("\n echo \"Primary Analysis Started\"\n")
	shellScriptFile.write("# Step 1: Generate VCFs\n\n")
	# Here I'm going to write out a Bam2gVCF for each of the bam files listed.  I'll just add a HC.g.vcf tag to the bam file names
	for i in range(1,len(BAMS)+1,1):
		shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -nct $NSLOTS  -T HaplotypeCaller  \\\n")
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
		shellScriptFile.write("-R $GENOME_FASTA -I $WORKING_DIR$SAMPLE%d_BAM \\\n"%(i))
		shellScriptFile.write("-o $WORKING_DIR${SAMPLE%d_BAM}.HC.vcf \n"%i)
	#Here I'll merge the gVCFs generated above. So we can just refer to them within this code as:  $WORKING_DIR$SAMPLE%d_BAM.HC.g.vcf



# Merge the GVCFs, assuming you started with BAMs and these were generated in Step 1
<<<<<<< HEAD
def MergeGVCF_withBAMLIST(shellScriptFile,BAMS):
	shellScriptFile.write("\n# Merge gVCFs\n\n")
	shellScriptFile.write("$JAVA -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -T GenotypeGVCFs -nt $NSLOTS \\\n")
=======
def MergeGVCF_withBAMLIST():
	shellScriptFile.write("\n# Step 2: Merge gVCFs\n\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt $NSLOTS \\\n")
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	shellScriptFile.write("-R $GENOME_FASTA \\\n")
	for i in range(1,len(BAMS)+1,1):
		shellScriptFile.write("--variant $WORKING_DIR${SAMPLE%d_BAM}.HC.g.vcf \\\n"%i)
	shellScriptFile.write("-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
<<<<<<< HEAD

# Merge VCFs, assuming you started with BAMS and these were generated in Step 1
def MergeVCF_withBAMLIST(shellScriptFile,BAMS):
	shellScriptFile.write("\n# Merge VCFs\n\n")
	shellScriptFile.write("$JAVA -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -T CombineVariants \\\n")
=======
        shellScriptFile.write("\n#BGZIP that bad boy\n")
        shellScriptFile.write("#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
        shellScriptFile.write("#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz \n\n")

# Merge VCFs, assuming you started with BAMS and these were generated in Step 1
def MergeVCF_withBAMLIST():
	shellScriptFile.write("\n# Step 2: Merge VCFs\n\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T CombineVariants \\\n")
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	shellScriptFile.write("-R $GENOME_FASTA \\\n")
	for i in range(1,len(BAMS)+1,1):
		shellScriptFile.write("--variant $WORKING_DIR${SAMPLE%d_BAM}.HC.vcf \\\n"%i)
	shellScriptFile.write("-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
<<<<<<< HEAD

#If  you are starting with your VCF list, GVCFs
# This runs genotype GVCFs
def MergeGVCF_withVCFLIST(shellScriptFile,VCFS):
	shellScriptFile.write("\n# Merge gVCFs\n\n")
	shellScriptFile.write("$JAVA -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -T GenotypeGVCFs -nt 4 \\\n")
=======
	shellScriptFile.write("\n#BGZIP that bad boy\n")
        shellScriptFile.write("#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
        shellScriptFile.write("#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz \n\n")

#If  you are starting with your VCF list, GVCFs
# This runs genotype GVCFs
def MergeGVCF_withVCFLIST():
	shellScriptFile.write("\n# Step 2: Merge gVCFs\n\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 4 \\\n")
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	shellScriptFile.write("-R $GENOME_FASTA \\\n")
	for i in range(1,len(VCFS)+1,1):
		shellScriptFile.write("--variant $WORKING_DIR${SAMPLE%d_VCF} \\\n"%i)
	shellScriptFile.write("-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
<<<<<<< HEAD

# If you're starting with your VCF list, VCFs
# This runs combine variants
def MergeVCF_withVCFLIST(shellScriptFile,VCFS):
	shellScriptFile.write("\n# Merge VCFs\n\n")
	shellScriptFile.write("$JAVA -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -T CombineVariants \\\n")
=======
	shellScriptFile.write("\n#BGZIP that bad boy\n")
        shellScriptFile.write("#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
        shellScriptFile.write("#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz \n\n")

# If you're starting with your VCF list, VCFs
# This runs combine variants
def MergeVCF_withVCFLIST():
	shellScriptFile.write("\n# Step 2: Merge VCFs\n\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T CombineVariants \\\n")
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	shellScriptFile.write("-R $GENOME_FASTA \\\n")
	for i in range(1,len(VCFS)+1,1):
                shellScriptFile.write("--variant $WORKING_DIR${SAMPLE%d_VCF} \\\n"%i)
	shellScriptFile.write("-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
	shellScriptFile.write("\n#Get Rid of non-chr chromosomes\n")
<<<<<<< HEAD
=======
	shellScriptFile.write("\n#BGZIP that bad boy\n")
	shellScriptFile.write("#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf\n")
	shellScriptFile.write("#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz\n\n")
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	

# Doesn't matter if you have your own VCF list (either version) or a BAM list
# You'll still need to normalize your merged VCF
<<<<<<< HEAD
def MergedVCF2NormVCF(shellScriptFile):
	shellScriptFile.write("\n#  Normalize merged VCF, annotate with SNPeff\n\n")
	shellScriptFile.write("zless $VCF \\\n")
	shellScriptFile.write("\t| sed 's/ID=AD,Number=./ID=AD,Number=R/' \\\n")
	shellScriptFile.write("\t| $VT decompose -s - \\\n")
	shellScriptFile.write("\t| $VT normalize -r $GENOME_FASTA - \\\n")
	shellScriptFile.write("\t| $JAVA -Xmx10g -jar $SNPEFFJAR GRCh37.75 \\\n")
	shellScriptFile.write("\t| $BGZIP -c > $NORMVCF \n")
	shellScriptFile.write("$TABIX -p vcf $NORMVCF\n")

# Add a soft-filter to the merged, normalized, SNPeff VCF using BCFTools
def FilterVCF(shellScriptFile):
	shellScriptFile.write("\n# Filter Merged, normalized VCF\n\n")
	shellScriptFile.write("$BCFTOOLS filter \\\n")
	shellScriptFile.write("\t --include 'FORMAT/AD[*:1]>=5 && FORMAT/DP[*] < 600' \\\n")
=======
def MergedVCF2NormVCF():
	shellScriptFile.write("\n# Step 3: Normalize merged VCF\n\n")
	shellScriptFile.write("zless $VCF \\\n")
	shellScriptFile.write("\t| sed 's/ID=AD,Number=./ID=AD,Number=R/' \\\n")
	shellScriptFile.write("\t| /opt/tools/vt/vt decompose -s - \\\n")
	shellScriptFile.write("\t| /opt/tools/vt/vt normalize -r $GENOME_FASTA - \\\n")
	shellScriptFile.write("\t| java -Xmx10g -jar $SNPEFFJAR GRCh37.75 \\\n")
	shellScriptFile.write("\t| /opt/tools/tabix/bgzip -c > $NORMVCF \n")
	shellScriptFile.write("/opt/tools/tabix/tabix -p vcf $NORMVCF\n")



# Made modification here for filter
def FilterVCF():
	shellScriptFile.write("\n# Step 4: Filter Merged, normalized VCF\n\n")
	shellScriptFile.write("/opt/tools/bcftools-1.8/bin/bcftools filter \\\n")
	shellScriptFile.write("\t --include 'FORMAT/AD[*:1]>=5 && FORMAT/DP[*] < 600 && FORMAT/AD[' \\\n")
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	shellScriptFile.write("\t -m + \\\n")
	shellScriptFile.write("\t -s + \\\n")
	shellScriptFile.write("\t -O z \\\n")
	shellScriptFile.write("\t --output $NORMFILTERVCF \\\n")
	shellScriptFile.write("\t $NORMVCF \n\n")
<<<<<<< HEAD
	shellScriptFile.write("$TABIX $NORMFILTERVCF \\\n\n")
	
# This is basic, just running the VCFAnno command on your input merged VCF
def RunVCFAnno(shellScriptFile,DatabaseDir):
	shellScriptFile.write("\n# VCFAnno - Turn your VCF file into an annotated VCF file\n")
	shellScriptFile.write('$VCFANNO -lua $ANNOTVARDIR/VCFAnno/custom.lua \\\n')
	shellScriptFile.write('-p $NSLOTS \'-base-path\' %s \\\n'%DatabaseDir)
	shellScriptFile.write('$ANNOTVARDIR/VCFAnno/VCFAnno_Config_20190321_GAC.toml \\\n')
	shellScriptFile.write('$NORMFILTERVCF > $ANNOVCF \n\n')

# NOTE: If you want to add certain things as --a-ok make sure you add them here, otherwise they may error on the creation of the mysqlDB
def VCF2DB(shellScriptFile):
	shellScriptFile.write("\n# VCF2DB - Turn your annotated VCF file into a GEMINI DB\n\n")
	shellScriptFile.write('python $VCF2DB \\\n')
=======
	shellScriptFile.write("/opt/tools/tabix/tabix $NORMFILTERVCF \\\n\n")

	
# This is basic, just running the VCFAnno command on your input merged VCF
def RunVCFAnno():
	shellScriptFile.write("\n# Step 5: VCFAnno - Turn your VCF file into an annotated VCF file\n")
	shellScriptFile.write('/opt/tools/vcfanno/vcfanno -lua /mnt/causes-vnx1/PIPELINES/AnnotateVariants/VCFAnno/custom.lua \\\n')
	shellScriptFile.write('-p $NSLOTS \\\n')
	shellScriptFile.write('/mnt/causes-vnx1/PIPELINES/AnnotateVariants/VCFAnno/VCFAnno_Config_20190131_GAC.toml \\\n')
	shellScriptFile.write('$NORMFILTERVCF > $ANNOVCF \n\n')

# NOTE: If you want to add certain things as --a-ok make sure you add them here, otherwise they may error on the creation of the mysqlDB
def VCF2DB():
	shellScriptFile.write("\n# Step 6: VCF2DB - Turn your annotated VCF file into a GEMINI DB\n\n")
	shellScriptFile.write('python /opt/tools/vcf2db/vcf2db.py \\\n')
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	shellScriptFile.write('--expand gt_quals --expand gt_depths --expand gt_alt_depths --expand gt_ref_depths --expand gt_types \\\n')
	shellScriptFile.write(' --a-ok InHouseDB_AC  --a-ok in_segdup --a-ok AF --a-ok AC --a-ok AN --a-ok MLEAC --a-ok MLEAF --a-ok gnomad_genome_hom_global --a-ok gnomad_genome_hom_afr --a-ok gnomad_genome_hom_amr --a-ok gnomad_genome_hom_asj --a-ok gnomad_genome_hom_eas --a-ok gnomad_genome_hom_fin --a-ok gnomad_genome_hom_nfe --a-ok gnomad_genome_hom_oth --a-ok gnomad_exome_hom_global --a-ok gnomad_exome_hom_afr --a-ok gnomad_exome_hom_amr --a-ok gnomad_exome_hom_asj --a-ok gnomad_exome_hom_eas --a-ok gnomad_exome_hom_fin --a-ok gnomad_exome_hom_nfe --a-ok gnomad_exome_hom_oth --a-ok cpg_island --a-ok common_pathogenic --a-ok cse-hiseq --a-ok DS --a-ok ConfidentRegion \\\n')
	shellScriptFile.write('$ANNOVCF $PED_FILE $GEMINIDB \n')

<<<<<<< HEAD

def Singleton_RenameVCF2MergedVCF(shellScriptFile):
	shellScriptFile.write("\n# Step 1-2: rename your VCF to merged VCF\n")
	shellScriptFile.write("cp $WORKING_DIR$SAMPLE1_VCF $WORKING_DIR${FAMILY_ID}.merged.hc.vcf\n")


def Main():

	# Get the arguments 
	args=GetArgs()


	# Print out some sanity checks
	# Make sure you have the correct working Directory.
	print "The working directory is: %s"%args.workingDir
	print "Family you're working with: %s"%args.Family

=======
def Main():


	args=GetArgs()
	print args
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
#######################
### Write the Script###
#######################

	#This is where we'll begin writing the script, starting with the header information
<<<<<<< HEAD
	shellScriptFile = open('%s%s_Bam2Gemini.sh'%(args.workingDir,args.Family),'w')
	print "You are generating this script: %s%s_Bam2Gemini.sh"%(args.workingDir,args.Family)
=======
	shellScriptFile = open('%s%s_Bam2Gemini.sh'%(workingDir,args.Family),'w')
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
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
<<<<<<< HEAD
	shellScriptFile.write('#PBS -o %s%s.o\n'%(args.workingDir,args.Family))
	shellScriptFile.write('#PBS -e %s%s.e\n'%(args.workingDir,args.Family))
=======
	shellScriptFile.write('#PBS -o %s%s.o\n'%(workingDir,args.Family))
	shellScriptFile.write('#PBS -e %s%s.e\n'%(workingDir,args.Family))
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
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
	
<<<<<<< HEAD
	# Make some variables for tools which are used within this pipeline
	shellScriptFile.write("# Define Tool paths. If they are in your path, simply change these full filepaths to only be the final command\n")
	shellScriptFile.write("# For example: Change BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools to be BCFTOOLS=bcftools if it's in your path \n\n")
	shellScriptFile.write("ANNOTVARDIR=%s\n"%args.AnnotateVariantsDir)
	shellScriptFile.write("SNPEFFJAR=/opt/tools/snpEff/snpEff.jar\n")
	shellScriptFile.write("BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools\n")
	shellScriptFile.write("VCFANNO=/opt/tools/vcfanno/vcfanno\n")
	shellScriptFile.write("VCF2DB=/opt/tools/vcf2db/vcf2db.py\n")
	shellScriptFile.write("GATKJAR=/opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar\n")
	shellScriptFile.write("JAVA=/opt/tools/jdk1.7.0_79/bin/java\n")
	shellScriptFile.write("BGZIP=/opt/tools/tabix/bgzip\n")
	shellScriptFile.write("TABIX=/opt/tools/tabix/tabix\n")
	shellScriptFile.write("VT=/opt/tools/vt/vt\n\n")
		

	# Add variable names to the shell script
	shellScriptFile.write("# Define variables for what you are working with\n\n")

	#Set the variables for working directory, sampleID, and fastq full filepath
	shellScriptFile.write("FAMILY_ID=\'%s\'\n"%args.Family)
	shellScriptFile.write("WORKING_DIR=\'%s\'\n"%args.workingDir)
=======
	
	
		
	#Set the variables for working directory, sampleID, and fastq full filepath
	shellScriptFile.write("FAMILY_ID=\'%s\'\n"%args.Family)
	shellScriptFile.write("WORKING_DIR=\'%s\'\n"%workingDir)
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
	if args.GENOME=='hg19':
		shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa\'\n")
	elif args.GENOME=='GSC':
		shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa\'\n")
	else:
		print "You did not choose a viable genome version, choose either GSC or hg19"
		sys.exit()
	shellScriptFile.write("PED_FILE=$WORKING_DIR/%s\n"%args.PED)
	shellScriptFile.write("TMPDIR=${WORKING_DIR}tmpdir/\n")
	shellScriptFile.write("mkdir $TMPDIR\n")
<<<<<<< HEAD
	shellScriptFile.write("GEMINIDB=$WORKING_DIR${FAMILY_ID}.db\n")
	shellScriptFile.write("VCF=$WORKING_DIR${FAMILY_ID}.merged.hc.vcf\n")
	shellScriptFile.write("NORMVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcf.gz\n")
	shellScriptFile.write("NORMFILTERVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.filter.vcf.gz\n")
	shellScriptFile.write('ANNOVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcfanno.vcf.gz \n')

	
	#The BAM files from our list, if running in the BAM mode
	if args.BAMLIST:
		BAMS = args.BAMLIST.split(',')
		print "These are the BAMS you are working with"
		for i in range(1,len(BAMS)+1,1):
			shellScriptFile.write("SAMPLE%d_BAM=%s\n"%(i,BAMS[i-1]))
			print BAMS[i-1]
	# The VCF files from our list, if running in the VCF mode
	if args.VCFLIST:
		VCFS = args.VCFLIST.split(',')
		for i in range(1,len(VCFS)+1,1):
			shellScriptFile.write("SAMPLE%d_VCF=%s\n"%(i,VCFS[i-1]))
=======
	
	#The BAM files from our list, if running in the BAM mode
	if args.BAMLIST:
		for i in range(1,len(BAMS)+1,1):
			shellScriptFile.write("SAMPLE%d_BAM=%s\n"%(i,BAMS[i-1]))
	
	# The VCF files from our list, if running in the VCF mode
	if args.VCFLIST:
		for i in range(1,len(VCFS)+1,1):
			shellScriptFile.write("SAMPLE%d_VCF=%s\n"%(i,VCFS[i-1]))
		# Make sure you have the correct working Directory.
		workingDir = args.workingDir
		print "The working directory is: %s"%workingDir
		print "Family you're working with: %s"%args.Family
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
		
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
				sys.exit()
	
<<<<<<< HEAD
#### Call functions to populate the commands in the script

	# Run in either VCF or BAM mode


	if args.VCFLIST:
	# if running in GVCF mode, then you want to joint genotype from multiple GVCFs
		if args.vcftype == 'GVCF':
			if args.Singleton:
				print "You cannot use singleton with GVCF mode"
				sys.exit()
			print "You have chosen a script which will do the following:"
			print "1) Joint Genotype the GVCFs"
			MergeGVCF_withVCFLIST(shellScriptFile,VCFS)
			print "2) Normalize and run SNPeff"
			MergedVCF2NormVCF(shellScriptFile)
			print "3) Filter with BCFTools soft-filter"
			FilterVCF(shellScriptFile)
			print "4) Annotate with VCFAnno"
			RunVCFAnno(shellScriptFile,args.DatabaseDir)
			print "5) Create a gemini database"
			VCF2DB(shellScriptFile)
	# if running in the VCF mode, I assume you have separate VCFs whcih you want to merge into a single merged.vcf
		elif args.vcftype == 'VCF':
			print "You have chosen a script which will do the following:" 
			print "1) Merge the VCFs into a single merged VCF"
			if args.Singleton:
				Singleton_RenameVCF2MergedVCF(shellScriptFile)
			else:
				MergeVCF_withVCFLIST(shellScriptFile,VCFS)
			print "2) Normalize and run SNPeff"
			MergedVCF2NormVCF(shellScriptFile)
			print "3) Filter with BCFTools soft-filter"
			FilterVCF(shellScriptFile)
			print "4) Annotate with VCFAnno"
			RunVCFAnno(shellScriptFile,args.DatabaseDir)
			print "5) Create a gemini database"
			VCF2DB(shellScriptFile)
	
	elif args.BAMLIST:
		if args.vcftype == 'GVCF':
			if args.Singleton:
                                print "You cannot use singleton with GVCF mode"
                                sys.exit()
			print "You have chosen a script which will do the following:"
			print "1) Make GVCF files for each of your BAM inputs"
			Bam2GVCF(shellScriptFile,BAMS)
			print "2) Joint Genotype the GVCF files you created"
			MergeGVCF_withBAMLIST(shellScriptFile,BAMS)
			print "3) Normalize and run SNPeff"
			MergedVCF2NormVCF(shellScriptFile)
			print "4) Filter with BCFTools soft-filter"
			FilterVCF(shellScriptFile)
			print "5) Annotate with VCFAnno"
			RunVCFAnno(shellScriptFile,args.DatabaseDir)
			print "6) Create a gemini database"
			VCF2DB(shellScriptFile)
	
		elif args.vcftype == 'VCF':
			print "You have chosen a script which will do the following:"
			print "1) Call variants from BAM files to make VCFs"
			Bam2VCF(shellScriptFile,BAMS)
			print "2) Merge the VCFs into a single merged VCF"
			if args.Singleton:
				Singleton_RenameVCF2MergedVCF(shellScriptFile)
			else:
				MergeVCF_withVCFLIST(shellScriptFile,VCFS)
			print "3) Normalize and run SNPeff"
	                MergedVCF2NormVCF(shellScriptFile)
			print "4) Filter with BCFTools soft-filter"
	                FilterVCF(shellScriptFile)
			print "5) Annotate with VCFAnno"
			RunVCFAnno(shellScriptFile,args.DatabaseDir)
			print "6) Create a GEMINI database"
			VCF2DB(shellScriptFile)
	


if __name__=="__main__":
	Main()
=======
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
		if args.vcftype == 'GVCF':
			MergeGVCF_withVCFLIST()
			MergedVCF2NormVCF()
			FilterVCF()
			RunVCFAnno()
			VCF2DB()
		elif args.vcftype == 'VCF':
			MergeVCF_withVCFLIST()
			MergedVCF2NormVCF()
			FilterVCF()
			RunVCFAnno()
			VCF2DB()
	
	elif args.BAMLIST:
		if args.vcftype == 'GVCF':
			Bam2GVCF()
			MergeGVCF_withBAMLIST()
			MergedVCF2NormVCF()
			FilterVCF()
			RunVCFAnno()
			VCF2DB()
	
		elif args.vcftype == 'VCF':
			Bam2VCF()
			MergeGVCF_withBAMLIST()
	                MergedVCF2NormVCF()
	                FilterVCF()
			RunVCFAnno()
			VCF2DB()
	
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3
