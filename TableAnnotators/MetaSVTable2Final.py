#########################################
# Author: Phillip Richmond              #
# Contact: phillip.a.richmond@gmail.com #
# License: open source GNU              #
#########################################

# The purpose of this python script is to take in a gemini table file, and add some annotations. 
# The annotations are all gene-based at this point, and are simple to add from text files. The heavy lifting of finding the gene has already been done
# This script also assumes only 1 gene being in the annotation column. For cases of overlapping genes, or intergenic variants, this script will not handle those annotations well (as of April 23rd 2018)

###########
# Imports #
###########

import sys, argparse
import re


########################
# Function Definitions #
########################

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--InFile",help="This is a file that has come as an output from the GEMINI query/built-in functions (e.g. de novo)",type=str)
    parser.add_argument("-o","--OutFile",help="Output Gemini Table file with Additional Annotations",type=str)
    args = parser.parse_args()
    infilename=args.InFile
    outfilename=args.OutFile
    return infilename,outfilename


# This script relies on a few different files in order to add some manual annotations. Those are defined in the main function at the bottom.


## Get a dictionary for the RefSeq Gene Summaries, gene is key, value is summary
def GetSummaryDict(SummaryFileName):
    SummaryFile = open(SummaryFileName,'r')
    SummaryForGenes = {}
    for line in SummaryFile:
        Number,Symbol,Summary = line.strip('\n').split("\t")
        SummaryForGenes[Symbol]=Summary
    return SummaryForGenes

## Get a dictionary for the OMIM gene map, gene is key, value is Mim#
def GetOMIM_Gene2MimDict(OMIM_Gene2MimFileName):
    OMIMfile = open(OMIM_Gene2MimFileName,'r')
    GenesToMim = {}
    for line in OMIMfile:
        if line[0]=='#':
            continue
        cols=line.strip('\n').split('\t')
        mimNumber= cols[0]
        if cols[1]!='gene':
            continue
        GeneSymbol=cols[3]
        GenesToMim[GeneSymbol]=mimNumber
    return GenesToMim

## Get a dictionary for the OMIM gene to phenotypes, gene is key, phenotypes are the values, with multiple phenos split by '|'
def GetOMIM_Gene2PhenoDict(OMIM_Gene2PhenoFileName):
	OMIMfile = open(OMIM_Gene2PhenoFileName,'r')
	GeneToPheno = {}
	for line in OMIMfile:
		if line[0]=='#':
			continue
		cols = line.strip('\n').split('\t')
		gene = cols[0]
		phenos = "|".join(cols[1:])
		GeneToPheno[gene]=phenos	
	return GeneToPheno

#### 
# RVIS file looks like this:
#GENE	ALL_0.01%	%ALL_0.01%
#A1BG	-0.35	29.49

def GetRVISDict(RVISFILE):
	infile = open(RVISFILE,'r')
	Gene2RVIS = {}
	headerline = infile.readline()
	for line in infile:
		gene,rvis_score,rvis_pct = line.strip('\n').split('\t')
		Gene2RVIS[gene] = [rvis_score,rvis_pct]
	print Gene2RVIS
	return Gene2RVIS


# PLI file looks like this:
# gene	pLI	pRec
#AGRN	0.17335234048116	0.826647657926747

def GetPLIDict(PLIFILE):
	infile = open(PLIFILE,'r')
	Gene2PLI = {}
	headerline = infile.readline()
	for line in infile:
		gene,pli_score,pli_pct = line.strip('\n').split('\t')
		Gene2PLI[gene] = [pli_score,pli_pct]
	print Gene2PLI
	return Gene2PLI


# This function takes in the dictionaries desribed above, reads in the gemini infile, and outputs the gemini outfile
def AddColumnsToTable(GeminiInFileName,GeminiOutFileName,Gene2Pheno,Gene2Mim,GeneSummary,Gene2PLI,Gene2RVIS):
	infile = open(GeminiInFileName,'r')
	outfile = open(GeminiOutFileName,'w')
	header = infile.readline()
	header = header.strip('\n')
	outfile.write("%s\tOMIM_Entry\tOMIM_Phenotypes\tRVIS_Score\tRVIS_Pct\tpLI_Score\tpLI_Pct\tGeneSummary\n"%header)
	for line in infile:
		line = line.strip('\n')
		cols = line.split("\t")
		gene = cols[5]
		# Check for Omim phenotype, add if there, if not make it '.'
		if Gene2Pheno.has_key(gene):
			omim_pheno=Gene2Pheno[gene]
		else:
			omim_pheno='.'
		# Check for Omim gene, add hyperlink for excel if there, if not make it '.'
		if Gene2Mim.has_key(gene):
			omim_key=Gene2Mim[gene]
			omim_hyperlink='=HYPERLINK(\"http://www.omim.org/entry/%s\")'%omim_key
		else:
			omim_hyperlink='.'
	
		if Gene2PLI.has_key(gene):
			pLI_score = Gene2PLI[gene][0]
			pLI_pct = Gene2PLI[gene][1]
		else:
			pLI_score = '.'
			pLI_pct = '.'		
		
		if Gene2RVIS.has_key(gene):
			rvis_score = Gene2RVIS[gene][0]
			rvis_pct = Gene2RVIS[gene][1]
		else:
			rvis_score = '.'
			rvis_pct = '.'			

 		# Check for gene summary, if there, add it, if not, make it '.' 
		if GeneSummary.has_key(gene):
			gene_summary=GeneSummary[gene]
		else:
			gene_summary='.'

		outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(line,omim_pheno,omim_hyperlink,rvis_score,rvis_pct,pLI_score,pLI_pct,gene_summary))

		#reset variables
		gene_summary = '.'
		omim_pheno = '.'
		omim_hyperlink='.'
		rvis_score = '.'
		rvis_pct = '.'
		pLI_score = '.'
		pLI_pct = '.'			

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

# These are hard coded locations for database files. These files are small, so they are all just text files.
# They contain gene-based information
    SummaryFileName = '/mnt/causes-data01/data/Databases/RefSeqGene_Summaries_270316.txt'
    Gene2MimFileName = '/mnt/causes-data01/data/Databases/OMIM_mim2gene'
    Gene2Disease = '/mnt/causes-data01/data/Databases/OMIM_phenotype_genelist'
    PLI = '/mnt/causes-data01/data/Databases/TOLERANCE/PLI_March2016.txt'
    RVIS = '/mnt/causes-data01/data/Databases/TOLERANCE/RVIS_March2016.txt'


    # Read in the annotations
    GeneSummaries = GetSummaryDict(SummaryFileName)
    Gene2Mim = GetOMIM_Gene2MimDict(Gene2MimFileName)
    Gene2Pheno = GetOMIM_Gene2PhenoDict(Gene2Disease)
    Gene2PLI = GetPLIDict(PLI)
    Gene2RVIS = GetRVISDict(RVIS)

	
    GeminiInfile,GeminiOutfile=GetOptions()
    #GeminiInfile='/mnt/causes-data01/data/RICHMOND/AnnotateVariants/T008_compoundHet.txt'
    #GeminiOutfile='/mnt/causes-data01/data/RICHMOND/AnnotateVariants/T008_compoundHet_annotated.txt'
 
    AddColumnsToTable(GeminiInfile,GeminiOutfile,Gene2Pheno,Gene2Mim,GeneSummaries,Gene2PLI,Gene2RVIS)




    sys.exit()



