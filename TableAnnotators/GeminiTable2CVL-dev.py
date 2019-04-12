#########################################
# Author: Phillip Richmond              #
# Contact: phillip.a.richmond@gmail.com #
# License: open source GNU              #
#########################################

# The purpose of this python script is to take in a gemini table file, and add some annotations. 
# The annotations are all gene-based at this point, and are simple to add from text files. The heavy lifting of finding the gene has already been done
# This script also assumes only 1 gene being in the annotation column. For cases of overlapping genes, or intergenic variants, this script will not handle those annotations well (as of April 23rd 2018)
# March 21st 2019 - Added OE scores, added gnomad hyperlink (https://gnomad.broadinstitute.org/variant/20-2465304-T-A)

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
	#print Gene2RVIS
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
	#print Gene2PLI
	return Gene2PLI

def GetHPODict(HPOFILE):
	infile = open(HPOFILE,'r')
	Gene2HPO={}
	headerline = infile.readline()
	for line in infile:
		entrez,gene,pheno,hpo = line.strip('\n').split('\t')
		if Gene2HPO.has_key(gene):
			Gene2HPO[gene].append("%s|%s"%(pheno,hpo))
		else:
			Gene2HPO[gene]=["%s|%s"%(pheno,hpo)]
	#print Gene2HPO
	return Gene2HPO

def GetMESHOPDict(MESHOPFILE):
        infile = open(MESHOPFILE,'r')
        Gene2MESHOP={}
        for line in infile:
                gene,meshop = line.strip('\n').split('\t')
                if Gene2MESHOP.has_key(gene):
                        Gene2MESHOP[gene].append(meshop)
                else:
                        Gene2MESHOP[gene]=[meshop]
        #print Gene2MESHOP
        return Gene2MESHOP	


# Added 20190321
def MakeGnomADHyperlink(chrom,pos,ref,alt):
	base = 'https://gnomad.broadinstitute.org/variant/'
	if 'chr' in chrom:
		chrom = chrom[3:]
	gnomad_hyperlink='=HYPERLINK(\"%s%s-%s-%s-%s\")'%(base,chrom,pos,ref,alt)
	return gnomad_hyperlink


def MakeVarCardsHyperlink(chrom,pos,ref,alt):
	base = 'http://159.226.67.237/sun/varcards/search/variant/'
	if 'chr' in chrom:
                chrom = chrom[3:]
        varcards_hyperlink='=HYPERLINK(\"%s%s-%s-%s-%s\")'%(base,chrom,pos,ref,alt)
        return varcards_hyperlink

# OE score
def MakeOEDict(OEFILE):
	infile = open(OEFILE,'r')
	Gene2OE = {}
	headerline = infile.readline()
	for line in infile:
		cols = line.strip('\n').split('\t')
		gene = cols[0]
		oe_mis = cols[4]
		oe_lof = cols[23]
		pLI = cols[20]
		Gene2OE[gene] = [oe_mis,oe_lof,pLI]
	return Gene2OE

# FLAGS - Just a list of genes to 'flag'
def MakeFLAGSList(FLAGSFILE):
	infile = open(FLAGSFILE,'r')
	FLAGS_GeneList = []
	for line in infile:
		FLAGS_GeneList.append(line.strip('\n'))
	return FLAGS_GeneList

# Gene name mapping, this is a 1-to-1
def MakeGeneNameDict(GENENAMEFILE):
	infile = open(GENENAMEFILE,'r')
        Gene2NAME={}
        for line in infile:
                gene,name = line.strip('\n').split('\t')
                Gene2NAME[gene] = name
        return Gene2NAME


# Gene alias mapping, multiple lines may have multiple aliases, 
# Will store an array here for each gene, and collapse below (similar to HPO/MESHOP)
def MakeGeneAliasDict(GENEALIASFILE):
	infile = open(GENEALIASFILE,'r')
        Gene2ALIAS={}
        for line in infile:
                gene,alias = line.strip('\n').split('\t')
		if Gene2ALIAS.has_key(gene):
                	Gene2ALIAS[gene].append(alias)
		else:
			Gene2ALIAS[gene] = [alias]
        return Gene2ALIAS



# This function takes in the dictionaries desribed above, reads in the gemini infile, and outputs the gemini outfile
# NOTE: THIS IS HARD CODED TO A SPECIFIC TABLE FORMAT
def AddColumnsToTable(GeminiInFileName,GeminiOutFileName,Gene2Pheno,Gene2Mim,GeneSummary,Gene2PLI,Gene2RVIS,Gene2HPO,Gene2MESHOP,Gene2OE,FLAGS_GeneList,Gene2NAME,Gene2ALIAS,ArrayOrderFile):
	infile = open(GeminiInFileName,'r')
	outfile = open(GeminiOutFileName,'w')
	header = infile.readline()
	header = header.strip('\n')
	
	# newline,gnomad_hyperlink,varcards_hyperlink,omim_pheno,omim_hyperlink,flags,rvis_score,rvis_pct,pLI_score,oe_mis_score,oe_lof_score,hpo,meshop,gene_name,gene_alias,gene_summary
	outfile.write("%s\tgnomAD_Hyperlink\tVarCards_Hyperlink\tOMIM_Phenotypes\tOMIM_Entry\tFLAGS\tRVIS_Score\tRVIS_Pct\tpLI_Score\tOE_Missense\tOE_LoF\tHPO\tMeSHOP\tGene_Name\tGene_Alias\tGeneSummary\n"%header)
	for line in infile:
		line = line.strip('\n')
		cols = line.split("\t")
		gene = cols[5]
		chrom = cols[0]
		pos = cols[2]
		ref = cols[3]
		alt = cols[4]

		# Additions 20190321
		# add gnomad hyperlink (20190321)
		gnomad_hyperlink = MakeGnomADHyperlink(chrom,pos,ref,alt)

		# add varcards hyperlink 
		varcards_hyperlink = MakeVarCardsHyperlink(chrom,pos,ref,alt)

		# Add flags
		if gene in FLAGS_GeneList:
			flags = 'FLAGS Gene'
		else:
			flags = '.'

		# Gene Name
		if Gene2NAME.has_key(gene):
                        gene_name=Gene2NAME[gene]
                else:
                        gene_name='.'

		# Gene Alias
		if Gene2ALIAS.has_key(gene):
                        gene_alias=';'.join(Gene2ALIAS[gene])
                else:
                        gene_alias='.'

		# PLI, OE_missense, OE_LoF
		if Gene2OE.has_key(gene):
			pLI_score = Gene2OE[gene][2]
			oe_mis_score = Gene2OE[gene][0]
			oe_lof_score = Gene2OE[gene][1]
		else:
			pLI_score = '.'
			oe_mis_score = '.'
			oe_lof_score = '.'	
		
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

		# HPO		
		if Gene2HPO.has_key(gene):
			hpo=";".join(Gene2HPO[gene])
		else:
			hpo = '.'	
		
		# RVIS	
		if Gene2RVIS.has_key(gene):
			rvis_score = Gene2RVIS[gene][0]
			rvis_pct = Gene2RVIS[gene][1]
		else:
			rvis_score = '.'
			rvis_pct = '.'			

		# MESHOP
		if Gene2MESHOP.has_key(gene):
                        meshop=";".join(Gene2MESHOP[gene])
                else:
                        meshop = '.'	

 		# Check for gene summary, if there, add it, if not, make it '.' 
		if GeneSummary.has_key(gene):
			gene_summary=GeneSummary[gene]
		else:
			gene_summary='.'

		# Fix some formatting stuff for Excel
		# Exon number
		cols[6]="\'%s"%cols[6]
		# Gene name
		#cols[5]="\'%s"%cols[5]

		# join the cols
		newline = "\t".join(cols)
		
#		ReOrderCols(ArrayOrderFile,cols)
		
		outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(newline,gnomad_hyperlink,varcards_hyperlink,omim_pheno,omim_hyperlink,flags,rvis_score,rvis_pct,pLI_score,oe_mis_score,oe_lof_score,hpo,meshop,gene_name,gene_alias,gene_summary))

		#reset variables
		gene_summary = '.'
		gene_name = '.'
		gene_alias = '.'
		gnomad_hyperlink = '.'
		varcards_hyperlink = '.'
		omim_pheno = '.'
		omim_hyperlink='.'
		flags = '.'
		rvis_score = '.'
		rvis_pct = '.'
		pLI_score = '.'
		oe_mis_score = '.'
		oe_lof_score = '.'
		hpo = '.'
		meshop = '.'
		

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

# These are hard coded locations for database files. These files are small, so they are all just text files.
# They contain gene-based information
        SummaryFileName = '/mnt/causes-vnx1/DATABASES/RefSeqGene_Summaries_270316.txt'
        Gene2MimFileName = '/mnt/causes-vnx1/DATABASES/OMIM_mim2gene'
        Gene2Disease = '/mnt/causes-vnx1/DATABASES/OMIM_phenotype_genelist'
        PLI = '/mnt/causes-vnx1/DATABASES/TOLERANCE/PLI_March2016.txt'
        RVIS = '/mnt/causes-vnx1/DATABASES/TOLERANCE/RVIS_March2016.txt'
        MESHOP = '/mnt/causes-vnx1/DATABASES/gene2pubmedBG-hum-gene2pubmed-gene-mesh-p_ONLYMESHDISEAS_P-valuecorrected_withGeneSymbols.txt'
        HPO = '/mnt/causes-vnx1/DATABASES/ALL_SOURCES_FREQUENT_FEATURES_genes_to_phenotype.txt'
	ALIAS = '/mnt/causes-vnx1/DATABASES/GeneNameMapping/HGNC_approved_symbol_and_alias_symbol.txt'
	NAME = '/mnt/causes-vnx1/DATABASES/GeneNameMapping/HGNC_approved_symbol_and_approved_name.txt'
	OE = '/mnt/causes-vnx1/DATABASES/OE/gnomad.v2.1.1.lof_metrics.by_gene.txt'
	FLAGSFILE = '/mnt/causes-vnx1/DATABASES/FLAGS/FLAGS_genes__12920_2014_64_MOESM4_ESM.txt'
	ArrayOrderFile = '/mnt/causes-vnx1/PIPELINES/AnnotateVariants/TableAnnotators/TemplateHeaderOrder.txt'

    # Read in the annotations
        GeneSummaries = GetSummaryDict(SummaryFileName)
        Gene2Mim = GetOMIM_Gene2MimDict(Gene2MimFileName)
        Gene2Pheno = GetOMIM_Gene2PhenoDict(Gene2Disease)
        Gene2PLI = GetPLIDict(PLI)
        Gene2RVIS = GetRVISDict(RVIS)
        Gene2HPO = GetHPODict(HPO)
        Gene2MESHOP = GetMESHOPDict(MESHOP)
	Gene2ALIAS = MakeGeneAliasDict(ALIAS)
	FLAGS_GeneList =  MakeFLAGSList(FLAGSFILE)
	Gene2NAME = MakeGeneNameDict(NAME)
	Gene2OE = MakeOEDict(OE)
        GeminiInfile,GeminiOutfile=GetOptions()
 
        AddColumnsToTable(GeminiInfile,GeminiOutfile,Gene2Pheno,Gene2Mim,GeneSummaries,Gene2PLI,Gene2RVIS,Gene2HPO,Gene2MESHOP,Gene2OE,FLAGS_GeneList,Gene2NAME,Gene2ALIAS,ArrayOrderFile)

        sys.exit()


