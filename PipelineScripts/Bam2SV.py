import sys, os, argparse

#####################################
# Author:   Phillip Richmond        #
# Contact:  prichmond@cmmt.ubc.ca   #
# Open source GNU licensing         #
#####################################

<<<<<<< HEAD
########
#Bam2SV#
########

# The purpose of this script is to take in a BAM or a set of BAM files, and genotype them with structural variant calling tools
# This tool can work with exome or genome
=======
############
#Bam2Gemini#
############
>>>>>>> 25ea7d23b876b086e1e0cba8a19e94ff20d835c3


##################
### Initialize ###
##################

if len(sys.argv) < 2:
	print "Re-run with the -h option"
	print "Typical Running Command:"
	print "python Bam2Gemini.py  -d /mnt/causes-data04/PROCESS/GENOME_TIDEX/T014_new/ -p 16 -m 40G \\"
	print "-P T014.ped -F T014  \\"
	print "-B T014-1_BWAmem_dupremoved_realigned.sorted.bam,T014-2_BWAmem_dupremoved_realigned.sorted.bam,T014-3_BWAmem_dupremoved_realigned.sorted.bam\n"

	print "If running singleton:"
	print "python Bam2Gemini.py  -d /mnt/causes-vnx2/TIDE/PROCESS/EXOME_TIDEX/Metabolomics/TIDEX071/ -p 12 -m 40G \\"
        print "-P TIDEX071.ped -F TIDEX017 \\"
        print "-V TIDEX071_BWAmem_dupremoved_realigned.sorted.bam\n"

	sys.exit()



# Read in your arguments
parser = argparse.ArgumentParser()
parser.add_argument("-d","--workingDir",help="Working directory on your filesystem",required=True)
parser.add_argument("-F","--Family",help="The Family Identifier",required=True)
parser.add_argument("-p","--processors",help="Choose the number of processors for this job",type=int,required=True)
parser.add_argument("--cnv",help="Run CNV calling and annotation",action='store_true',default=False)
parser.add_argument("--sv",help="Run SV calling and annotation",action='store_true',default=False)
parser.add_argument("-m","--memory",help="Choose the memory needed for this job",required=True)
parser.add_argument("-P","--PED",help="The PED File for this family",required=True)
parser.add_argument("-B","--BAMLIST",help="Comma separated list of BAMS to be processed.  These BAMs should correspond to the identifiers inside the Ped File")
parser.add_argument("-E","--Email",help="Email address",type=str)
parser.add_argument("-G","--GENOME",help="Which Genome version do you want to use? Options are GSC || hg19",required=True)
#parser.add_argument("-C","--Config",help="Config file which points at locations for tool executables",required=True)
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
else:
	print "You must supply with BAM files"
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
	shellScriptFile.write("GENOME_FASTA=\'/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa\'\n")
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
else:
	print "You must provide BAM files"
	sys.exit()


#########################################################################################

# CNVnator - CNV Calling
def CNVNATOR():
        shellScriptFile.write("#Running CNVNATOR Windowsize 100\n\n")
        shellScriptFile.write("#Defining Variables\n")
        shellScriptFile.write("BAM=${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
        shellScriptFile.write("WIN=100\n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR${BAM}.root -genome $GENOME_FASTA -tree $WORKING_DIR$BAM \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -his $WIN -d $CHROM \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -stat $WIN \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -partition $WIN \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -call $WIN > $WORKING_DIR/${SAMPLE_ID}_CNVnatorCall_$WIN \n")

        shellScriptFile.write("#Running CNVNATOR Windowsize 1000\n\n")
        shellScriptFile.write("WIN=1000\n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -his $WIN -d $CHROM \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -stat $WIN \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -partition $WIN \n")
        shellScriptFile.write("/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -call $WIN > $WORKING_DIR/${SAMPLE_ID}_CNVnatorCall_$WIN \n")


# ERDS - CNV Calling
def ERDS():
        shellScriptFile.write("\n# ERDS \n\n")
        shellScriptFile.write("perl /opt/tools/erds1.1/erds_pipeline.pl \\\n")
        shellScriptFile.write("\t-b $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \\\n")
        shellScriptFile.write("\t-v $WORKING_DIR${SAMPLE_ID}_dupremoved_realigned_HaplotypeCaller.vcf \\\n")
        shellScriptFile.write("\t-o $WORKING_DIR${SAMPLE_ID}_ERDS/ \\\n")
        shellScriptFile.write("\t-r $GENOME_FASTA \\\n")




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


# Run in either VCF or BAM mode
if args.BAMLIST:
	if args.CNV:
		ERDS()
		CNVNATOR()	
		if args.annotate:
			return
	if args.SV:
		SMOOVE()
		if args.annotate:
			return

	sys.exit()

		
