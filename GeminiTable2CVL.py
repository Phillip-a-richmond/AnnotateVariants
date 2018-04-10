#########################################
# Author: Phillip Richmond              #
# Contact: phillip.a.richmond@gmail.com #
# License: open source GNU              #
#########################################

# The purpose of this python script is to take in a gemini table file, and add some annotations. 

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


# This script relies on a few different files in order to add some manual annotations, including:
## RefSeq Gene Summaries
def GetSummaryDict(SummaryFileName):
    SummaryFile = open(SummaryFileName,'r')
    SummaryForGenes = {}
    for line in SummaryFile:
        Number,Symbol,Summary = line.strip('\n').split("\t")
        SummaryForGenes[Symbol]=Summary
    return SummaryForGenes

## OMIM Gene Map
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

def AddColumnsToTable(GeminiInFileName,GeminiOutFileName,Gene2Pheno,Gene2Mim,GeneSummary):
	infile = open(GeminiInFileName,'r')
	outfile = open(GeminiOutFileName,'w')
	header = infile.readline()
	header = header.strip('\n')
	outfile.write("%s\tOMIM_Entry\tOMIM_Phenotypes\tGeneSummary\n"%header)
	for line in infile:
		line = line.strip('\n')
		cols = line.split("\t")
		gene = cols[5]
		if Gene2Pheno.has_key(gene):
			omim_pheno=Gene2Pheno[gene]
		else:
			omim_pheno='.'
		if Gene2Mim.has_key(gene):
			omim_key=Gene2Mim[gene]
			omim_hyperlink='=HYPERLINK(\"http://www.omim.org/entry/%s\")'%omim_key
		else:
			ommim_hyperlink='.'
		if GeneSummary.has_key(gene):
			gene_summary=GeneSummary[gene]
		else:
			gene_summary='.'
		
		
		outfile.write("%s\t%s\t%s\t%s\n"%(line,omim_pheno,omim_hyperlink,gene_summary))

		#reset variables
		gene_summary = '.'
		omim_pheno = '.'
		omim_hyperlinke='.'

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    SummaryFileName = '/mnt/causes-data01/data/Databases/RefSeqGene_Summaries_270316.txt'
    Gene2MimFileName = '/mnt/causes-data01/data/Databases/OMIM_mim2gene'
    Gene2Disease = '/mnt/causes-data01/data/Databases/OMIM_phenotype_genelist'
    
    # Read in the annotations
    GeneSummaries = GetSummaryDict(SummaryFileName)
    Gene2Mim = GetOMIM_Gene2MimDict(Gene2MimFileName)
    Gene2Pheno = GetOMIM_Gene2PhenoDict(Gene2Disease)
	
    GeminiInfile,GeminiOutfile=GetOptions()
    #GeminiInfile='/mnt/causes-data01/data/RICHMOND/AnnotateVariants/T008_compoundHet.txt'
    #GeminiOutfile='/mnt/causes-data01/data/RICHMOND/AnnotateVariants/T008_compoundHet_annotated.txt'
 
    AddColumnsToTable(GeminiInfile,GeminiOutfile,Gene2Pheno,Gene2Mim,GeneSummaries)
    sys.exit()



