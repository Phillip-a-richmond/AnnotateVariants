import sys, os, argparse

#####################################
# Author:   Phillip Richmond        #
# Contact:  prichmond@cmmt.ubc.ca   #
# Open source GNU licensing         #
#####################################

############
#Bam2Gemini#
############


##################
### Initialize ###
##################

if len(sys.argv) < 2:
	print "Re-run with the -h option"
	print "Typical Running Command (From Bams, GVCF mode):"
	print "python Bam2Gemini.py  -d /mnt/causes-data04/PROCESS/GENOME_TIDEX/T014_new/ -p 16 -m 40G \\"
	print "-P T014.ped -F T014 -v GVCF \\"
	print "-B T014-1_BWAmem_dupremoved_realigned.sorted.bam,T014-2_BWAmem_dupremoved_realigned.sorted.bam,T014-3_BWAmem_dupremoved_realigned.sorted.bam\n"


	print "Typical Running Command (From GVCFs):"
        print "python Bam2Gemini.py  -d /mnt/causes-data04/PROCESS/GENOME_TIDEX/T014_new/ -p 16 -m 40G \\"
        print "-P T014.ped -F T014 -v GVCF \\"
        print "-V T014-1_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf,T014-2_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf,T014-3_BWAmem_dupremoved_realigned_HaplotypeCaller.g.vcf\n"

	print "Typical Running Command, duo (From VCFs):"
        print "python Bam2Gemini.py  -d /mnt/causes-data04/PROCESS/GENOME_TIDEX/T274_new/ -p 16 -m 40G \\"
        print "-P T274.ped -F T274 -v VCF \\"
        print "-V tidex1006_BWAmem_dupremoved_realigned_HaplotypeCaller.vcf,tidex1007_BWAmem_dupremoved_realigned_HaplotypeCaller.vcf\n"

	sys.exit()



# Read in your arguments
parser = argparse.ArgumentParser()
parser.add_argument("-d","--workingDir",help="Working directory on your filesystem",required=True)
parser.add_argument("-F","--Family",help="The Family Identifier",required=True)
parser.add_argument("-p","--processors",help="Choose the number of processors for this job",type=int,required=True)
parser.add_argument("-m","--memory",help="Choose the memory needed for this job",required=True)
parser.add_argument("-P","--PED",help="The PED File for this family",required=True)
parser.add_argument("-B","--BAMLIST",help="Comma separated list of BAMS to be processed.  These BAMs should correspond to the identifiers inside the Ped File")
parser.add_argument("-E","--Email",help="Email address",type=str)
parser.add_argument("-v","--vcftype",help="The type of VCF.  If these are GVCFs that need to be merged, or separately called VCFs. REQUIRED option.",required=True)
parser.add_argument("-V","--VCFLIST",help="Comma separated list of VCF files.  Set either GVCF or VCF with -v to know merging option.")
parser.add_argument("-G","--GENOME",help="Which Genome version do you want to use? Options are GSC || hg19",required=True)
#parser.add_argument("-C","--Config",help="Config file which points at locations for tool executables",required=True)
#parser.add_argument("-Q","--QueryScript",help="Query script template",required=True)
args = parser.parse_args()




# Make sure you have the correct working Directory.
workingDir = args.workingDir
print "The working directory is: %s"%workingDir
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
		sys.exit()



# These are parameters for the scheduler
numProcessors = args.processors
Memory = args.memory
runTime = "240:00:00"


#######################
### Write the Script###
#######################

#This is where we'll begin writing the script, starting with the header information
shellScriptFile = open('%s%s_Bam2Gemini.sh'%(workingDir,args.Family),'w')
shellScriptFile.write('#!/bin/bash\n')

#Job name
shellScriptFile.write('#PBS -N %s_Bam2Gemini\n'%args.Family)
#Export environment variables
shellScriptFile.write('#PBS -V\n')
#set location for log files
shellScriptFile.write('#PBS -o %s%s.o\n'%(workingDir,args.Family))
shellScriptFile.write('#PBS -e %s%s.e\n'%(workingDir,args.Family))
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
shellScriptFile.write("WORKING_DIR=\'%s\'\n"%workingDir)
if args.GENOME=='hg19':
	sys.exit("hg19 genome location unclear")
	shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-data01/data/GENOMES/hg19/FASTA/hg19.fa\'\n")
elif args.GENOME=='GSC':
	shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa\'\n")
else:
	print "You did not choose a viable genome version, choose either GSC or hg19"
	sys.exit()
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



#########################################################################################

# Generate GVCFs from BAM files
def Bam2GVCF():
	shellScriptFile.write("\n echo \"Primary Analysis Started\"\n")
	shellScriptFile.write("# Step 1: Generate gVCFs\n\n")
	# Here I'm going to write out a Bam2gVCF for each of the bam files listed.  I'll just add a HC.g.vcf tag to the bam file names
	for i in range(1,len(BAMS)+1,1):
	        shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar \\\n")
	        shellScriptFile.write(" -T HaplotypeCaller -nct 4 --emitRefConfidence GVCF \\\n")
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
def Bam2VCF():
	shellScriptFile.write("\n echo \"Primary Analysis Started\"\n")
	shellScriptFile.write("# Step 1: Generate VCFs\n\n")
	# Here I'm going to write out a Bam2gVCF for each of the bam files listed.  I'll just add a HC.g.vcf tag to the bam file names
	for i in range(1,len(BAMS)+1,1):
		shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -nct $NSLOTS  -T HaplotypeCaller  \\\n")
		shellScriptFile.write("-R $GENOME_FASTA -I $WORKING_DIR$SAMPLE%d_BAM \\\n"%(i))
		shellScriptFile.write("-o $WORKING_DIR${SAMPLE%d_BAM}.HC.vcf \n"%i)
	#Here I'll merge the gVCFs generated above. So we can just refer to them within this code as:  $WORKING_DIR$SAMPLE%d_BAM.HC.g.vcf



# Merge the GVCFs, assuming you started with BAMs and these were generated in Step 1
def MergeGVCF_withBAMLIST():
	shellScriptFile.write("\n# Step 2: Merge gVCFs\n\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt $NSLOTS \\\n")
	shellScriptFile.write("-R $GENOME_FASTA \\\n")
	for i in range(1,len(BAMS)+1,1):
		shellScriptFile.write("--variant $WORKING_DIR${SAMPLE%d_BAM}.HC.g.vcf \\\n"%i)
	shellScriptFile.write("-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
        shellScriptFile.write("\n#BGZIP that bad boy\n")
        shellScriptFile.write("#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
        shellScriptFile.write("#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz \n\n")

# Merge VCFs, assuming you started with BAMS and these were generated in Step 1
def MergeVCF_withBAMLIST():
	shellScriptFile.write("\n# Step 2: Merge VCFs\n\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T CombineVariants \\\n")
	shellScriptFile.write("-R $GENOME_FASTA \\\n")
	for i in range(1,len(BAMS)+1,1):
		shellScriptFile.write("--variant $WORKING_DIR${SAMPLE%d_BAM}.HC.vcf \\\n"%i)
	shellScriptFile.write("-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
	shellScriptFile.write("\n#BGZIP that bad boy\n")
        shellScriptFile.write("#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
        shellScriptFile.write("#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz \n\n")

#If  you are starting with your VCF list, GVCFs
# This runs genotype GVCFs
def MergeGVCF_withVCFLIST():
	shellScriptFile.write("\n# Step 2: Merge gVCFs\n\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 4 \\\n")
	shellScriptFile.write("-R $GENOME_FASTA \\\n")
	for i in range(1,len(VCFS)+1,1):
		shellScriptFile.write("--variant $WORKING_DIR${SAMPLE%d_VCF} \\\n"%i)
	shellScriptFile.write("-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
	shellScriptFile.write("\n#BGZIP that bad boy\n")
        shellScriptFile.write("#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
        shellScriptFile.write("#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz \n\n")

# If you're starting with your VCF list, VCFs
# This runs combine variants
def MergeVCF_withVCFLIST():
	shellScriptFile.write("\n# Step 2: Merge VCFs\n\n")
	shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -Djava.io.tmpdir=$TMPDIR -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T CombineVariants \\\n")
	shellScriptFile.write("-R $GENOME_FASTA \\\n")
	for i in range(1,len(VCFS)+1,1):
                shellScriptFile.write("--variant $WORKING_DIR${SAMPLE%d_VCF} \\\n"%i)
	shellScriptFile.write("-o $WORKING_DIR${FAMILY_ID}.merged.hc.vcf \n")
	shellScriptFile.write("\n#Get Rid of non-chr chromosomes\n")
	shellScriptFile.write("\n#BGZIP that bad boy\n")
	shellScriptFile.write("#/opt/tools/tabix/bgzip $WORKING_DIR${FAMILY_ID}.merged.hc.vcf\n")
	shellScriptFile.write("#/opt/tools/tabix/tabix $WORKING_DIR${FAMILY_ID}.merged.hc.vcf.gz\n\n")
	

# Doesn't matter if you have your own VCF list (either version) or a BAM list
# You'll still need to normalize your merged VCF
def MergedVCF2NormVCF():
	shellScriptFile.write("\n# Step 3: Normalize merged VCF\n\n")
	shellScriptFile.write("# Define some variables\n\n")
	shellScriptFile.write("SNPEFFJAR=/opt/tools/snpEff/snpEff.jar\n")
	shellScriptFile.write("GEMINIDB=$WORKING_DIR${FAMILY_ID}.db\n")
	shellScriptFile.write("VCF=$WORKING_DIR${FAMILY_ID}.merged.hc.vcf\n")
	shellScriptFile.write("NORMVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcf.gz\n")
	shellScriptFile.write("zless $VCF \\\n")
	shellScriptFile.write("\t| sed 's/ID=AD,Number=./ID=AD,Number=R/' \\\n")
	shellScriptFile.write("\t| /opt/tools/vt/vt decompose -s - \\\n")
	shellScriptFile.write("\t| /opt/tools/vt/vt normalize -r $GENOME_FASTA - \\\n")
	shellScriptFile.write("\t| java -Xmx10g -jar $SNPEFFJAR GRCh37.75 \\\n")
	shellScriptFile.write("\t| /opt/tools/tabix/bgzip -c > $NORMVCF \n")
	shellScriptFile.write("/opt/tools/tabix/tabix -p vcf $NORMVCF\n")

def FilterVCF():
	shellScriptFile.write("\n# Step 4: Filter Merged, normalized VCF\n\n")
	shellScriptFile.write("NORMFILTERVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.filter.vcf.gz\n")
	shellScriptFile.write("/opt/tools/bcftools-1.8/bin/bcftools filter \\\n")
	shellScriptFile.write("\t --include 'FORMAT/AD[*:1]>=10 && FORMAT/DP[*] < 300' \\\n")
	shellScriptFile.write("\t -m + \\\n")
	shellScriptFile.write("\t -s + \\\n")
	shellScriptFile.write("\t -O z \\\n")
	shellScriptFile.write("\t --output $NORMFILTERVCF \\\n")
	shellScriptFile.write("\t $NORMVCF \n\n")
	shellScriptFile.write("/opt/tools/tabix/tabix $NORMFILTERVCF \\\n\n")


def RunVCFAnno():
	shellScriptFile.write("\n# Step 5: VCFAnno - Turn your VCF file into an annotated VCF file\n")
	shellScriptFile.write('ANNOVCF=$WORKING_DIR${FAMILY_ID}.merged.hc.norm.vcfanno.vcf.gz \n')
	shellScriptFile.write('/opt/tools/vcfanno/vcfanno -lua /mnt/causes-vnx1/PIPELINES/AnnotateVariants/VCFAnno/custom.lua \\\n')
	shellScriptFile.write('-p $NSLOTS \\\n')
	shellScriptFile.write('/mnt/causes-vnx1/PIPELINES/AnnotateVariants/VCFAnno/VCFANNO_Config_PlusGNOMAD_PlusInHouse_SplitByPop_gnomAD_Exome_VNX.toml \\\n')
	shellScriptFile.write('$NORMFILTERVCF > $ANNOVCF \n\n')

# NOTE: If you want to add certain things as --a-ok make sure you add them here, otherwise they may error on the creation of the mysqlDB
def VCF2DB():
	shellScriptFile.write("\n# Step 6: VCF2DB - Turn your annotated VCF file into a GEMINI DB\n\n")
	shellScriptFile.write('python /opt/tools/vcf2db/vcf2db.py \\\n')
	shellScriptFile.write('--expand gt_quals --expand gt_depths --expand gt_alt_depths --expand gt_ref_depths --expand gt_types \\\n')
	shellScriptFile.write(' --a-ok InHouseDB_AC  --a-ok in_segdup --a-ok AF --a-ok AC --a-ok AN --a-ok MLEAC --a-ok MLEAF --a-ok gnomad_genome_hom_global --a-ok gnomad_genome_hom_afr --a-ok gnomad_genome_hom_amr --a-ok gnomad_genome_hom_asj --a-ok gnomad_genome_hom_eas --a-ok gnomad_genome_hom_fin --a-ok gnomad_genome_hom_nfe --a-ok gnomad_genome_hom_oth --a-ok gnomad_exome_hom_global --a-ok gnomad_exome_hom_afr --a-ok gnomad_exome_hom_amr --a-ok gnomad_exome_hom_asj --a-ok gnomad_exome_hom_eas --a-ok gnomad_exome_hom_fin --a-ok gnomad_exome_hom_nfe --a-ok gnomad_exome_hom_oth --a-ok cpg_island --a-ok common_pathogenic --a-ok cse-hiseq --a-ok DS --a-ok ConfidentRegion \\\n')
	shellScriptFile.write('$ANNOVCF $PED_FILE $GEMINIDB \n')

# This function will add the gemini build command from a merged haplo vcf
# This is deprecated now 
def MergedNormVCF2GeminiDB():
	shellScriptFile.write("\n# Step 4: Merged, Normalized VCF 2 Gemini DB\n\n")
	shellScriptFile.write("/opt/tools/gemini/bin/gemini load --cores $NSLOTS \\\n")
	shellScriptFile.write("-t snpEff -v $NORMFILTERVCF -p $PED_FILE --tempdir $TMPDIR $GEMINIDB \n")
# This is deprecated now
def AddGNOMAD2GeminiDB():
	shellScriptFile.write("\n# Step 5: Add gnomAD\n\n")
	shellScriptFile.write("genome_VCF=/mnt/causes-vnx1/DATABASES/GNOMAD/gnomad.genomes.r2.0.2.sites.wholeGenome.norm.vcf.gz\n")
	shellScriptFile.write("/opt/tools/gemini/bin/gemini annotate -f $genome_VCF -a extract -e AF,Hom -c aaf_gnomAD_genome_all,gnomAD_genome_num_hom_alt -t float,integer -o max,max $GEMINIDB\n")
	shellScriptFile.write('echo "UPDATE variants SET aaf_gnomAD_genome_all = -1.0, gnomAD_genome_num_hom_alt = -1 where aaf_gnomAD_genome_all is NULL;" | sqlite3 $GEMINIDB\n')


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

