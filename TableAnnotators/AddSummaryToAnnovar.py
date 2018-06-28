import sys


if len(sys.argv) < 3:
	print "Usage: python AddSummaryToAnnovar.py <in.annovar> <out.summaryadded.annovar>"
	sys.exit()


infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')

#The refseq ID stuff
SummaryFile = open('/mnt/causes-data01/data/Databases/RefSeqGene_Summaries_270316.txt','r')
SummaryForGenes = {}
for line in SummaryFile:
	Number,Symbol,Summary = line.strip('\n').split("\t")
	SummaryForGenes[Symbol]=Summary

#The OMIM Stuff
OMIMfile = open('/mnt/causes-data01/data/Databases/OMIM_mim2gene','r')
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
#print GenesToMim

header = infile.readline().strip('\n')
outfile.write("%s\tSummary\tOMIMhyperlink\n"%header)
for line in infile:
	cols = line.strip('\n').split('\t')
	Genes = cols[6]
	MimNumber = "."
	Summary = "."
	for each in Genes.split(','):
		if SummaryForGenes.has_key(each):
			Summary = SummaryForGenes[each]
		if GenesToMim.has_key(each):
			MimNumber = GenesToMim[each]
	cols.append(Summary)
	if len(MimNumber) > 2:
		hyperlink = "=HYPERLINK(\"http://www.omim.org/entry/%s\")"%MimNumber
	else:
		hyperlink=MimNumber
	cols.append(hyperlink)
	newline="\t".join(cols)
	outfile.write("%s\n"%newline)

	




