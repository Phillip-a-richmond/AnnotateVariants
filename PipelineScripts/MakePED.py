import sys, os, argparse

#####################################
# Author:   Phillip Richmond        #
# Contact:  prichmond@cmmt.ubc.ca   #
# Open source GNU licensing         #
#####################################

##########
#Make PED#
##########

# This tool will make a ped file from input arguments

##################
### Initialize ###
##################

if len(sys.argv) < 2:
	print "Re-run with the -h option"
	print "Typical Running Command:"
	print "python MakePED.py --proband TIDEX1000_BWAmem,male,affected --father TIDEX1002_BWAmem,male,unaffected --mother TIDEX1001_BWAmem,female,unaffected --family T272 -O T272.ped"
	print "WARNING: Only works for single-family set ups, and is currently limited to one sibling, and does not use ethnicity information."
	sys.exit()


# Parse options with Argparse
parser = argparse.ArgumentParser()
parser.add_argument("-p","--proband",help="Proband info in the form of: SampleID,male|female,affected|unaffected",required=True)
parser.add_argument("-s","--sibling",help="Sibling info in the form of: SampleID,male|female,affected|unaffected")
parser.add_argument("-f","--father",help="Father info in the form of: SampleID,male|female,affected|unaffected")
parser.add_argument("-m","--mother",help="Mother info in the form of: SampleID,male|female,affected|unaffected")
parser.add_argument("-O","--outfile",help="Output PED file to write to",required=True)
parser.add_argument("-F","--FamilyID",help="Family identifier",required=True)
args = parser.parse_args()


# Write out PED file 
outfile = open(args.outfile,'w')

# A really crude translation dictionary, just to make life a bit easier
TranslateDict = {'affected':2,'unaffected':1,'male':1,'female':2}

# Write the header
outfile.write("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tethnicity\n")



# Add the father
if args.father:
	fatherCols = args.father.split(',')
	fatherID = fatherCols[0]
	outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(args.FamilyID,fatherCols[0],'-9','-9',TranslateDict[fatherCols[1]],TranslateDict[fatherCols[2]],'-9'))
else: 
	fatherID = '-9'

# Add the mother
if args.mother:
	motherCols = args.mother.split(',')
	motherID = motherCols[0]
	outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(args.FamilyID,motherCols[0],'-9','-9',TranslateDict[motherCols[1]],TranslateDict[motherCols[2]],'-9'))
else:
	motherID = '-9'


# Add the proband
probandCols = args.proband.split(',')
outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(args.FamilyID,probandCols[0],fatherID,motherID,TranslateDict[probandCols[1]],TranslateDict[probandCols[2]],'-9'))

# Write out sibling
if args.sibling:
	siblingCols = args.sibling.split(',')
	outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(args.FamilyID,siblingCols[0],fatherID,motherID,TranslateDict[siblingCols[1]],TranslateDict[siblingCols[2]],'-9'))


