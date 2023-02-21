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
    parser.add_argument("-i","--InFile",help="This is a file that has come as an output from the GEMINI query/built-in functions (e.g. de novo)",type=str,required=True)
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
        Symbol,Number,Summary = line.strip('\n').split("\t")
        SummaryForGenes[Symbol]=Summary
    #print(SummaryForGenes)
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
                cols = line.strip('\n').split('\t')
                entrez = cols[0]
                gene = cols[1]
                hpo = cols[2]
                pheno = cols[3]
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
        end = '?dataset=gnomad_r3'
	if 'chr' in chrom:
		chrom = chrom[3:]
	gnomad_hyperlink='=HYPERLINK(\"%s%s-%s-%s-%s%s\")'%(base,chrom,pos,ref,alt,end)
	return gnomad_hyperlink

def MakeClinvarHyperlink(clinvarID):
        base = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/'
        clinvar_hyperlink = '=HYPERLINK(\"%s%s\")'%(base,clinvarID)
        return clinvar_hyperlink

def MakeGeneCardsHyperlink(gene):
	base = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='
        genecards_hyperlink='=HYPERLINK(\"%s%s\")'%(base,gene)
        return genecards_hyperlink

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


# The point of this function is to re-order the columns according to some order file
# Essentially, I want to read in this header column array, find the value of each index, and then take an array of values and pull the correct indices 
# repopulate a new array
def ReOrderCols(ORDERFILE,INARRAYHEADER,INARRAYVALS):


	infile = open(ORDERFILE,'r')
	# This is the array of final ordering, based on the input ordering file
	FinalArrayOrder = []
	
	for line in infile:
		FinalArrayOrder.append(line.strip('\n'))
	print FinalArrayOrder

	# Next, I'll read in the array header to a dictionary, where a key is e.g. 'gene' or 'chrom'
	# and the value is the position within the existing array
	HeaderDict = {}
	for i in range(len(INARRAYHEADER)):
		HeaderDict[INARRAYHEADER[i]] = i

	print HeaderDict

	# I'm going to build an empty array and repopulate it according to the final array order
	dummyarray = range(len(FinalArrayOrder))
	
	# Now I'll take the final array order, and for each item I'll get the index from the column
	for i in range(len(FinalArrayOrder)):
		val = FinalArrayOrder[i]
		dummyarray[i] = INARRAYVALS[HeaderDict[val]]
	print dummyarray
	ReorderedCols=dummyarray
	return FinalArrayOrder,ReorderedCols


# This function takes in the dictionaries desribed above, reads in the gemini infile, and outputs the gemini outfile
# NOTE: THIS IS HARD CODED TO A SPECIFIC TABLE FORMAT
def AddColumnsToTable(GeminiInFileName,GeminiOutFileName,Gene2Pheno,Gene2Mim,GeneSummary,Gene2PLI,Gene2RVIS,Gene2HPO,Gene2MESHOP,Gene2OE,FLAGS_GeneList,Gene2NAME,Gene2ALIAS,ArrayOrderFile):
	infile = open(GeminiInFileName,'r')
	outfile = open(GeminiOutFileName,'w')
	header = infile.readline()
	header = header.strip('\n')
	
	# Automated check to make sure there is something in this file, otherwise write to outfile: No Variants to Report 
	if len(header) < 1:
		outfile.write("No Variants to Report\n")
		return

	# This is the existing header, and I'll add the new features to this array
	headercols = header.split('\t')
	headercols.append('gnomad_hyperlink')
	headercols.append('genecards_hyperlink')
        headercols.append('clinvar_hyperlink')
	headercols.append('OMIM_Phenotypes')
	headercols.append('OMIM_Entry')
	headercols.append('FLAGS')
	headercols.append('RVIS_Score')
	headercols.append('RVIS_Pct')
	headercols.append('pLI_Score')
	headercols.append('OE_Missense')
	headercols.append('OE_LoF')
	headercols.append('HPO')
	headercols.append('MeSHOP')
	headercols.append('Gene_Name')
	headercols.append('Gene_Alias')
	headercols.append('GeneSummary')

	#print headercols
	#print headercols.index('gnomAD_Hyperlink')
	
	# a dummy test here
	cols = headercols
	FinalOrderHeader,ReorderedCols = ReOrderCols(ArrayOrderFile,headercols,cols)
	outfile.write("%s\n"%'\t'.join(FinalOrderHeader))
	#outfile.write("%s\tgnomAD_Hyperlink\tVarCards_Hyperlink\tOMIM_Phenotypes\tOMIM_Entry\tFLAGS\tRVIS_Score\tRVIS_Pct\tpLI_Score\tOE_Missense\tOE_LoF\tHPO\tMeSHOP\tGene_Name\tGene_Alias\tGeneSummary\n"%header)
	for line in infile:
		line = line.strip('\n')
		cols = line.split("\t")
		gene = cols[5]
		chrom = cols[0]
		pos = cols[2]
		ref = cols[3]
		alt = cols[4]
                # added 2021-05-10, just after gene to make it easy
                clinvarID = cols[6]
                # If there is no clinvarID, then this gets input as a NoneType, but need it to be a blank string
                if clinvarID == None:
                    clinvarID = '.'
                    cols[6] = '.'

		# Additions 20190321
		# add gnomad hyperlink (20190321)
		gnomad_hyperlink = MakeGnomADHyperlink(chrom,pos,ref,alt)

		# add genecards hyperlink  (2021-05-10)
		genecards_hyperlink = MakeGeneCardsHyperlink(gene)

                # add clinvar hyperlink (2021-05-10)
                clinvar_hyperlink = MakeClinvarHyperlink(clinvarID)

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
                        #print(gene_summary)
			gene_summary=GeneSummary[gene]
		else:
			gene_summary='.'

		# Fix some formatting stuff for Excel
		# Exon number
		cols[8]="\'%s"%cols[8]
		# Gene name
		cols[5]="\'%s"%cols[5]

		# join the cols
		newline = "\t".join(cols)
		
		cols.append(gnomad_hyperlink)
		cols.append(genecards_hyperlink)
                cols.append(clinvar_hyperlink)
		cols.append(omim_pheno)
		cols.append(omim_hyperlink)
		cols.append(flags)
		cols.append(rvis_score)
		cols.append(rvis_pct)
		cols.append(pLI_score)
		cols.append(oe_mis_score)
		cols.append(oe_lof_score)
		cols.append(hpo)
		cols.append(meshop)
		cols.append(gene_name)
		cols.append(gene_alias)
		cols.append(gene_summary)
		
		FinalOrderHeader,ReorderedCols = ReOrderCols(ArrayOrderFile,headercols,cols)
                # make empty cells '.'
                for i in range(len(ReorderedCols)):
                    if (ReorderedCols[i]=='') or (ReorderedCols[i]==' ') or (ReorderedCols[i]==None) or (ReorderedCols[i]=='None') or (ReorderedCols[i]=="'") or (ReorderedCols[i]=='.'):
                        ReorderedCols[i]='None'
		outfile.write("%s\n"%'\t'.join(ReorderedCols))

		#reset variables
		gene_summary = '.'
		gene_name = '.'
		gene_alias = '.'
		gnomad_hyperlink = '.'
		genecards_hyperlink = '.'
                clinvar_hyperlink = '.'
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
        SummaryFileName = '/mnt/common/DATABASES/GENERIC/GeneNameMapping/GeneSymbolGeneNameGeneSummary.txt'
        Gene2MimFileName = '/mnt/common/DATABASES/GENERIC/OMIM/OMIM_mim2gene'
        Gene2Disease = '/mnt/common/DATABASES/GENERIC/OMIM/OMIM_phenotype_genelist'
        PLI = '/mnt/common/DATABASES/GENERIC/PLI/PLI_March2016.txt'
        RVIS = '/mnt/common/DATABASES/GENERIC/PLI/RVIS_March2016.txt'
        MESHOP = '/mnt/common/DATABASES/GENERIC/MESHOPS/gene2pubmedBG-hum-gene2pubmed-gene-mesh-p_ONLYMESHDISEAS_P-valuecorrected_withGeneSymbols.txt'
        HPO = '/mnt/common/DATABASES/GENERIC/HPO/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt'
	ALIAS = '/mnt/common/DATABASES/GENERIC/GeneNameMapping/HGNC_approved_symbol_and_alias_symbol.txt'
	NAME = '/mnt/common/DATABASES/GENERIC/GeneNameMapping/HGNC_approved_symbol_and_approved_name.txt'
	OE = '/mnt/common/DATABASES/GENERIC/OE/gnomad.v2.1.1.lof_metrics.by_gene.txt'
	FLAGSFILE = '/mnt/common/DATABASES/GENERIC/FLAGS/FLAGS_genes__12920_2014_64_MOESM4_ESM.txt'
	ArrayOrderFile = '/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/TableAnnotators/TemplateHeaderOrder_2021-05-10.txt'

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



