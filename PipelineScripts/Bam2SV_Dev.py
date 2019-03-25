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
        parser.add_argument("--sv",help="Run SV calling and annotation",action='store_true',default=False)
        parser.add_argument("--mei",help="Run MEI (mobile element insertion) calling",action='store_true',default=False)
        parser.add_argument("--STR",help="Run STR (short tandem repeat) calling with STRetch",action='store_true',default=False)
	args = parser.parse_args()
	return args	


#########################################################################################
# SMOOVE SV Calling and Annotation with AnnotSV
def SMOOVE(shellScriptFile,BAMS):
	shellScriptFile.write("#Running SMOOVE\n\n")
	shellScriptFile.write("docker run -v /mnt:/mnt brentp/smoove smoove call \\\n")
	shellScriptFile.write("--name $SAMPLE_ID -p $NSLOTS -x \\\n")
	shellScriptFile.write("--outdir $WORKING_DIR/Smoove \\\n")
	shellScriptFile.write("--genotype --duphold -f $GENOME_FASTA \\\n")
	for i in range(1,len(BAMS)+1,1):
                shellScriptFile.write("${SAMPLE%d_BAM} \\\n"%i)

# CNV Calling with CNVnator
def CNVNATOR(shellScriptFile,BAMS):
        shellScriptFile.write("#Running CNVNATOR Windowsize 100\n\n")
        shellScriptFile.write("#Defining Variables\n")
        shellScriptFile.write("WIN=100\n")
	for i in range(1,len(BAMS)+1,1):
		shellScriptFile.write("BAM=${SAMPLE%d_BAM}\n"%i)
        	shellScriptFile.write("$CNVNATOR -root $WORKING_DIR${BAM}.root -genome $GENOME_FASTA -tree $WORKING_DIR$BAM \n")
        	shellScriptFile.write("$CNVNATOR -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -his $WIN -d $CHROM \n")
        	shellScriptFile.write("$CNVNATOR -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -stat $WIN \n")
        	shellScriptFile.write("$CNVNATOR -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -partition $WIN \n")
        	shellScriptFile.write("$CNVNATOR -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -call $WIN > $WORKING_DIR/${SAMPLE_ID}_CNVnatorCall_$WIN \n")

        shellScriptFile.write("#Running CNVNATOR Windowsize 1000\n\n")
        shellScriptFile.write("WIN=1000\n")
	for i in range(1,len(BAMS)+1,1):
                shellScriptFile.write("BAM=${SAMPLE%d_BAM}\n"%i)
       		shellScriptFile.write("$CNVNATOR -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -his $WIN -d $CHROM \n")
       		shellScriptFile.write("$CNVNATOR -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -stat $WIN \n")
       		shellScriptFile.write("$CNVNATOR -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -partition $WIN \n")
       		shellScriptFile.write("$CNVNATOR -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -call $WIN > $WORKING_DIR/${SAMPLE_ID}_CNVnatorCall_$WIN \n")

#CNV Calling with ERDS
def ERDS(shellScriptFile,BAMS,VCFS):
        shellScriptFile.write("\n# ERDS \n\n")
	for i in range(1,len(BAMS)+1,1):
                shellScriptFile.write("BAM=${SAMPLE%d_BAM}\n"%i)
	        shellScriptFile.write("perl $ERDS  \\\n")
	        shellScriptFile.write("\t-b $BAM \\\n")
	        shellScriptFile.write("\t-v $VCF \\\n")
	        shellScriptFile.write("\t-o $WORKING_DIR${BAM}_ERDS/ \\\n")
	        shellScriptFile.write("\t-r $GENOME_FASTA \\\n")

# SV Calling with Pindel
def Pindel(shellScriptFile,args):
        shellScriptFile.write("# Running Pindel \n\n")
        shellScriptFile.write("#Defining Variables\n")
	for i in range(1,len(BAMS)+1,1):
                shellScriptFile.write("BAM=${SAMPLE%d_BAM}\n"%i)
	        shellScriptFile.write("\n##Generate Empirical insert size stats\n")
	        shellScriptFile.write("$SAMTOOLS view -r $SAMPLE_ID $WORKING_DIR$BAM | tail -n+1000000 | python /opt/tools/lumpy/pairend_distro.py -r %s -X 4 -N 1000000 -o $WORKING_DIR${BAM}.histo > $WORKING_DIR${BAM}.insertStats\n"%args.readlength)
	        shellScriptFile.write("MEAN=`cat $WORKING_DIR${BAM}.insertStats | sed -E 's/\s+/,/' | cut -d, -f1 | sed -E 's/mean://' | xargs printf \"%.0f\"`\n")
	        shellScriptFile.write("python /mnt/causes-vnx1/PIPELINES/AnnotateVariants/pindel_config.py $WORKING_DIR$BAM $MEAN $WORKING_DIR$SAMPLE_ID \n")
	        shellScriptFile.write("/opt/tools/pindel-0.2.5b6/pindel --number_of_threads $NSLOTS \\\n")
	        shellScriptFile.write("-f $GENOME_FASTA \\\n")
	        shellScriptFile.write("-i $WORKING_DIR${BAM}_config.txt \\\n")
	        shellScriptFile.write("-c ALL -o $WORKING_DIR$SAMPLE_ID \\\n")
	        shellScriptFile.write("-M 4 \\\n")
	        shellScriptFile.write("-N -x 3 \\\n")
	        shellScriptFile.write("-I false -t false -r false \\\n")
	        shellScriptFile.write("-J /mnt/causes-vnx1/DATABASES/hg19_centromeres_telomeres.bed  \n\n")
	

# Mobile element insertion calling
def MEI(shellScriptFile):
        shellScriptFile.write("\n# Mobile Element Insertions\n")
        shellScriptFile.write("## Define Variables\n")
        shellScriptFile.write("ANALYSIS_DIR=${WORKING_DIR}MEI\n")
        shellScriptFile.write("MELT_DIR=/opt/tools/MELVTv2.1.5/\n")
        shellScriptFile.write("MEI_LIST=${MELT_DIR}/me_refs/1KGP_Hg19/mei_list.txt\n")
        shellScriptFile.write("GENE_ANNO=/opt/tools/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed\n\n")
        for i in range(1,len(BAMS)+1,1):
                shellScriptFile.write("BAM=${SAMPLE%d_BAM}\n"%i)
	shellScriptFile.write("# MELT Singleton \n\n")
        shellScriptFile.write("java -jar ${MELT_DIR}MELT.jar Single \\\n")
        shellScriptFile.write("-a -b hs37d5/NC007605 -c 8 -h $GENOME_FASTA \\\n")
        shellScriptFile.write("-bamfile $WORKING_DIR$$BAM \\\n")
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

###### STR calling ######
def STRetch_GRCh37(shellScriptFile):
        shellScriptFile.write("\n# Short Tandem Repeats\n")
	for i in range(1,len(BAMS)+1,1):
                shellScriptFile.write("BAM=${SAMPLE%d_BAM}\n"%i)
	        shellScriptFile.write("/mnt/causes-vnx1/PIPELINES/STRetch/tools/bin/bpipe run \\\n")
	        shellScriptFile.write("-p input_regions=/mnt/causes-vnx1/PIPELINES/STRetch/reference-data/GRCh37.simpleRepeat_period1-6_dedup.sorted.bed \\\n")
	        shellScriptFile.write("/mnt/causes-vnx1/PIPELINES/STRetch/pipelines/GRCh37_STRetch_wgs_bam_pipeline.groovy \\\n")
	        shellScriptFile.write("$WORKING_DIR${SAMPLE_ID}_dupremoved_realigned.sorted.bam \n")
	
def STRetch_hg19(shellScriptFile):
        shellScriptFile.write("\n# Short Tandem Repeats\n")
	for i in range(1,len(BAMS)+1,1):
                shellScriptFile.write("BAM=${SAMPLE%d_BAM}\n"%i)
	        shellScriptFile.write("/mnt/causes-vnx1/PIPELINES/STRetch/tools/bin/bpipe run \\\n")
	        shellScriptFile.write("-p input_regions=/mnt/causes-vnx1/PIPELINES/STRetch/reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed \\\n")
	        shellScriptFile.write("/mnt/causes-vnx1/PIPELINES/STRetch/pipelines/hg19_STRetch_wgs_bam_pipeline.groovy \\\n")
	        shellScriptFile.write("$WORKING_DIR$BAM \n")
	

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
	
	
	
		
	# Add variable names to the shell script
	shellScriptFile.write("# Define some variables\n\n")
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
                BAMS = args.BAMLIST.split(',')
                print "These are the BAMS you are working with"
                for i in range(1,len(BAMS)+1,1):
                        shellScriptFile.write("SAMPLE%d_BAM=%s\n"%(i,BAMS[i-1]))
                        print BAMS[i-1]
		
	# NOTE: BAM files are referred to within this pipeline simply as: SAMPLE#_BAM, where # is an integer for the BAM number


	

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
        shellScriptFile.write("CNVNATOR=/opt/tools/CNVnator/src/cnvnator\n")
        shellScriptFile.write("FASTQC=/opt/tools/FastQC/fastqc\n")
        shellScriptFile.write("ERDS=/opt/tools/erds1.1/erds_pipeline.pl\n")



#### Call functions to populate the commands in the script

	
	if args.BAMLIST:
		if (args.mei):
                        MEI(shellScriptFile)
                if args.cnv:
               		shellScriptFile.write("\necho \"CNV Analysis Started\"\n")
                        shellScriptFile.write("date\n")
                        CNVNATOR(shellScriptFile)
                        ERDS(shellScriptFile)
                if args.STR:
                        shellScriptFile.write("\necho \"STR calling\"\n")
                        if args.GENOME=='hg19':
                                STRetch_hg19(shellScriptFile)
                        elif args.GENOME=='GSC':
                                STRetch_GRCh37(shellScriptFile)
                if args.sv:
			
			
		
	else:
		print "You did not supply the BAMs you fool!"
		sys.exit()



