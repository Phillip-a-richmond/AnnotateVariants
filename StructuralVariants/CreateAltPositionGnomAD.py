import sys, argparse

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-I","--Infile",help="Input gnomAD SV bed file", type=str, required=True)
    parser.add_argument("-O","--OutPrefix",help="Output prefix for gnomAD altered files", type=str, required=True)
    args = parser.parse_args()
    infilename=args.Infile
    outfileprefix=args.OutPrefix
    return infilename,outfileprefix


# Purpose of this function is to take in the full gnomadSV bed file, extract the insertion calls, and then list the alternate coordinates for those insertion calls
# Reasoning behind this is that when we run MELT, it looks like the annotation is pinging the alternative loci, and not providing these useful annotation flags, e.g. an, ac, af
 
def ParseInputWrite_INS_shifted(infilename,outfileprefix):
    infile = open(infilename,'r')
    outfile = open("%s.ins.shifted.bed"%outfileprefix, 'w')
    # Make an outfile with chr* for the chrom names
    outfile2 = open("%s.ins.shifted.chr.bed"%outfileprefix, 'w')

    # Read through the file, grab the columns you want
    # Select the INS calls
    # write output to new file

    # Columns we want here are for the alternate mapping locus of the insertion calls
    # which in this case are shifted by...what seems to be the size of an LTR?
    # In our input, these columns are relevant:
    # shifted coords: chr2 (col[7]), pos2 (col[16]), end2 (col[10])
    # variant identifier: name (col[3])
    # SV length: svlen (col[30])
    # SV type: svtype (col[31])
    # allele freq: AN (col[34]), AC (col[35]), AF (col[36]), NHOMALT (col[40]), POPMAXAF (col[69])
    # So together we'll make a new file that has: 
    # chr2,pos2,end2,name,svlen,svtype,AN,AC,AF,NHOMALT,POPMAXAF
    for line in infile.readlines():
        cols = line.strip('\n').split('\t')
        chr2=cols[7]
        pos2=cols[16]
        end2=cols[10]
        name=cols[3]
        svlen=cols[30]
        svtype=cols[31]
        an=cols[35]
        ac=cols[36]
        af=cols[37]
        nhomalt=cols[41]
        popmaxaf=cols[70]
        # catch header
        if line[0] == '#':
            #print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(chr2,pos2,end2,name,svlen,svtype,an,ac,af,nhomalt,popmaxaf))
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chr2,pos2,end2,name,svlen,svtype,an,ac,af,nhomalt,popmaxaf))
        # keep only the insertions
        if svtype=='INS':
            #print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(chr2,pos2,end2,name,svlen,svtype,an,ac,af,nhomalt,popmaxaf))
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chr2,pos2,end2,name,svlen,svtype,an,ac,af,nhomalt,popmaxaf))
            outfile2.write('chr%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chr2,pos2,end2,name,svlen,svtype,an,ac,af,nhomalt,popmaxaf))



# Purpose of this function is to input the gnomadSV bed file, extract informative columns, and print the output.
# Reasoning behind this is that the annotations provided by AnnotSV don't include useful info like allele frequency, etc.
def ParseInputWrite_ALL_reformat(infilename,outfileprefix):
    infile = open(infilename,'r')
    outfile = open("%s.reformatted.bed"%outfileprefix,'w')
    outfile2 = open("%s.reformatted.chr.bed"%outfileprefix, 'w')


    # Read through the file, grab the columns you want
    # write output to new file

    # Columns we want here are primary position columns, sv len, and allele freqs
    # In our input, these columns are relevant:
    # coords: chr (col[0]), start (col[1]), end (col[2])
    # variant identifier: name (col[3])
    # SV length: svlen (col[30])
    # SV type: svtype (col[31])
    # allele freq: AN (col[34]), AC (col[35]), AF (col[36]), NHOMALT (col[40]), POPMAXAF (col[69])
    # So together we'll make a new file that has: 
    # chr2,pos2,end2,name,svlen,svtype,AN,AC,AF,NHOMALT,POPMAXAF
    for line in infile.readlines():
        cols = line.strip('\n').split('\t')
        chrom=cols[0]
        start=cols[1]
        end=cols[2]
        name=cols[3]
        svlen=cols[30]
        svtype=cols[31]
        an=cols[35]
        ac=cols[36]
        af=cols[37]
        nhomalt=cols[41]
        popmaxaf=cols[70]
        #print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(chrom,start,end,name,svlen,svtype,an,ac,af,nhomalt,popmaxaf))
        outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chrom,start,end,name,svlen,svtype,an,ac,af,nhomalt,popmaxaf))
        outfile2.write('chr%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chrom,start,end,name,svlen,svtype,an,ac,af,nhomalt,popmaxaf))



def Main():
    Infilename,Outfileprefix = GetOptions()
    print("Infile: %s\nOutfilePrefix: %s\n"%(Infilename,Outfileprefix))
    ParseInputWrite_INS_shifted(Infilename,Outfileprefix)
    ParseInputWrite_ALL_reformat(Infilename,Outfileprefix)

if __name__=="__main__":
    Main()



