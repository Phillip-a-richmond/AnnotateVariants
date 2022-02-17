#########################################
# Author: Phillip Richmond              #
# Contact: phillip.a.richmond@gmail.com #
# License: open source GNU              #
#########################################

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
    parser.add_argument("-i","--InFile",help="Input CVL",type=str,required=True)
    parser.add_argument("-o","--OutFile",help="Output CVL With GeneList filter", type=str, required=True)
    parser.add_argument("-g","--GeneList",help="List of genes, one gene per line", type=str, required=True)
    args = parser.parse_args()
    infilename=args.InFile
    outfilename=args.OutFile
    genelistfilename = args.GeneList
    return infilename,outfilename,genelistfilename


def ReadGeneListToArray(genelistfilename):
    genelistfile = open(genelistfilename,'r')
    GeneList = []
    for line in genelistfile:
        gene = line.strip('\n')
        GeneList.append(gene)
    return GeneList

def ParseCVLForGeneList(GeneList, infilename, outfilename):
    infile = open(infilename,'r')
    outfile = open(outfilename,'w')

    for line in infile:
        cols = line.strip('\n').split('\t')
        # keep header info with these 2 if statements
        if cols[0] == 'Nothing':
            outfile.write(line)
            continue
        if len(cols) < 10:
            outfile.write(line)
            continue
        # Now the only lines left should be those with actual genes in them, trimming the leading '
        gene = cols[1][1:]
        if gene in GeneList:
            print("Found a hit")
            print(gene)
            outfile.write(line)
        else:
            continue
        
		

#-------------#
# Main        #
#-------------#

def Main():
    InFileName,OutFileName,GeneListFileName = GetOptions()
    GeneList = ReadGeneListToArray(GeneListFileName)
    ParseCVLForGeneList(GeneList,InFileName,OutFileName)



if __name__ == "__main__":
    Main()
    sys.exit()



