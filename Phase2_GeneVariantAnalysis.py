# The purpose of this tool is to take the gene line, extract the variants, and make a table of them

#########################################
# Author: Phillip Richmond              #
# Contact: phillip.a.richmond@gmail.com #
# License: open source GNU              #
#########################################


# Each gene line looks like this:d
# GeneName  GeneMim#    chromosome  start   stop    cytolocation    computed_cyto   Gene Symbols    Gene Descriptive Name   Entrez Gene ID  Ensembl Gene ID Notes   Phenotype/Disease   MGI IMSR Info [List of available mice]  VariantList (Tab separated, chr_pos_ref_alt) [GRCh37/hg19]

# Example:
# PAX6  607108  chr11   31784791    31817960    11p13   11p13   PAX6, AN2, MGDA, FVH1, ASGD5    Paired box homeotic gene-6  5080    ENSG00000007372 mutation identified in 1 patient each with MDGA, COLBN, or COLB Aniridia, 106210 (3), Autosomal dominant; Anterior segment dysgenesis 5, multiple subtypes, 604229 (3); Cataract with late-onset corneal dystrophy, 106210 (3), Autosomal dominant; ?Coloboma of optic nerve, 120430 (3), Autosomal dominant; ?Coloboma, ocular, 120200 (3), Autosomal dominant; Foveal hypoplasia 1, 136520 (3), Autosomal dominant; Keratitis, 148190 (3), Autosomal dominant; ?Morning glory disc anomaly, 120430 (3), Autosomal dominant; Optic nerve hypoplasia, 165550 (3), Autosomal dominant    Pax6 (MGI:97490)    ['CMMR:0229 - 10-3-2_Pax6', 'EM:00008_Pax6', 'EM:00437_Pax6', 'EM:00440_Pax6', 'EM:00426_Pax6', 'EM:01331_Pax6', 'EM:01329_Pax6', 'EM:01328_Pax6', 'EM:01327_Pax6', 'EM:01326_Pax6', 'EM:01325_Pax6', 'EM:01324_Pax6', 'EM:01323_Pax6', 'EM:01322_Pax6', 'EM:01330_Pax6', 'HAR:479_Pax6', 'HAR:121_Pax6', 'HAR:1052_Pax6', 'HAR:231_Pax6', 'HAR:1236_Pax6', 'HAR:870_Pax6', 'HAR:1077_Pax6', 'HAR:1238_Pax6', 'HAR:907_Pax6', 'HAR:1589_Pax6', 'JAX:024578_Pax6', 'JAX:028032_Pax6', 'JAX:024688_Pax6', 'JAX:024688_Pax6', 'MMRRC:036955_Pax6', 'RBRC-GSC0045_Pax6', 'RBRC-GSC0031_Pax6', 'RBRC-GSC0014_Pax6', 'RBRC05945_Pax6', 'RBRC06513_Pax6', 'RBRC06457_Pax6', 'EM:00998_Pax6', 'EM:00011_Pax6', 'EM:01687_Pax6', 'EM:00440_Pax6', 'EM:00426_Pax6', 'EM:01686_Pax6', 'EM:01891_Pax6', 'HAR:479_Pax6', 'HAR:121_Pax6', 'HAR:1052_Pax6', 'HAR:231_Pax6', 'HAR:1236_Pax6', 'HAR:870_Pax6', 'HAR:1077_Pax6', 'HAR:1238_Pax6', 'HAR:907_Pax6', 'HAR:1589_Pax6', 'JAX:027971_Pax6', 'JAX:027971_Pax6', 'JAX:000391_Pax6', 'MMRRC:036955_Pax6', 'RBRC01273_Pax6']    chr11_31811483_T_A  chr11_31812270_T_C  chr11_31812383_G_C  chr11_31815069_G_A  chr11_31815335_G_A  chr11_31815343_A_G  chr11_31815345_C_T  chr11_31815620_C_G  chr11_31815627_G_A  chr11_31816238_G_A  chr11_31816247_G_A  chr11_31816253_G_A  chr11_31822356_G_A  chr11_31822380_G_A  chr11_31822385_A_T  chr11_31823109_G_T  chr11_31823159_G_A  chr11_31823264_G_A  chr11_31823267_T_A  chr11_31823267_T_G  chr11_31823275_C_A  chr11_31823314_C_A  chr11_31823316_G_A  chr11_31823316_G_T  chr11_31823441_A_T  chr11_31824281_G_A  chr11_31824317_G_C

###########
# Imports #
###########

import sys, argparse
import Conservation as conservationWithAA
import pybedtools
import re


########################
# Function Definitions #
########################


def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-G","--GeneFile",help="This is a file that has one gene line in it, copied and pasted from the online gene table.",type=str)
    parser.add_argument("-o","--outfile",help="Output variant table file",type=str)
    parser.add_argument("-W","--WindowSize",help="Window centered on variant",type=int,default=300)
    args = parser.parse_args()

    infile = open(args.GeneFile,'r')
    outfile = open(args.outfile,'w')
    windowsize=args.WindowSize
    return infile,outfile,windowsize


#given a nt sequence, return reverse compliment
def reverse_compliment(seq):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(seq) 
    letters = [basecomplement[base] for base in letters] 
    comp = ''.join(letters)
    revcomp = comp[::-1]
    return revcomp

#takes a motif, translates to regex, returns list with indices of positions within the string
def FindMotifInSeq(seq,motif):
    # define IUPAC dict
    IUPAC_dict = {'N':'[ACGT]','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','A':'A','T':'T','C':'C','G':'G'}
    # Build pattern using for loop for regex search
    pattern = ""
    for nt in motif:
        pattern = pattern+IUPAC_dict[nt]

    # Make sure we have all uppercase
    seq = seq.upper()
    revcompseq = reverse_compliment(seq) 

    #search in the forward orientation
    forwardhits = re.finditer(pattern,seq)
    forwardpositions =  [m.start(0) for m in forwardhits]

    #search in the reverse orientation
    reversehits = re.finditer(pattern,revcompseq)
    reversepositions = [m.start(0) for m in reversehits]
    
    return forwardpositions,reversepositions


# Parses the clinvar summary file, and produces 2 dictionaries that can be used to annotate variants in a downstream function
# Currently, the file that it parses, "variant_summary.txt" which I get from:
# ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
# has some entries in it that don't list the reference and alternate alleles on the variant line...
# so for those variants, I'm just using the chromosome and position to pull them out
# it's not great, but it works for now
def ReadClinvarVariantSummary(clinvarfilename):
    clinvarfile = open(clinvarfilename,'r')

    
    #dictionary where variants are stored using the rs tag as the key
    VARINFORS = {}
    
    #dictionary where I store information re: variants, stored with key as: chrom_pos_ref_alt
    VARINFO = {}

    # similar dictionary, but only for use when the ref-alt columns in the variant_summary.txt file are fucking wrong.  What kind of asshole does that?
    VARINFOSHORT = {}


    # read header line
    headerline = clinvarfile.readline()
    for line in clinvarfile:
        line=line.strip('\n')
        # skip the GRCh38 lines
        if "GRCh38" in line:
            continue
        if "Pathogenic" not in line:
            continue
        cols = line.split('\t')
        chromosome = 'chr%s'%cols[18]
        pos = cols[19]
        ref = cols[21]
        alt = cols[22]
        rsnum = "rs%s"%cols[9]
        name = cols[2]
        phenos = cols[12]
        phenos2 =cols[13]
        IDshort = "%s_%s"%(chromosome,pos)
        ID = "%s_%s_%s_%s"%(chromosome,pos,ref,alt)
        if len(ref) + len(alt) == 2:
            VARINFO[ID]=[rsnum,name,phenos,phenos2]
            continue
        # some lines don't have correct ref and alt listings where they should be in variant_summary
        VARINFOSHORT[IDshort]=[rsnum,name,phenos,phenos2]
        VARINFORS[rsnum]=[rsnum,name,phenos,phenos2]
    return VARINFO,VARINFOSHORT,VARINFORS



# The main program, that is taking in the variant info (From clinvar variant_summary.txt), the infile with the gene line, and a window size, and creates the table of variants
# The table will look something like this:
#chromosome,start,stop,position,ref,alt,seq,NT_conservation,AA_conservation,gnomAD_url,rsnum,name,phenos,phenos2,hitset,numhits,numCloseHits

def MakeVarDict(INFILE,OUTFILE,WindowSize,VarInfo,VarInfoShort,VarInfoRS,DataDir):

    # before we get into the fun bit, lets print out our header
    OUTFILE.write("Chromosome\tWindowStart\tWindowStop\tPosition\tReferenceAllele\tAlternateAllele\tSequence\tNT_Human\tAA_Human\t")
    OUTFILE.write("NT_Mouse_Seq\tNT_Mouse_Conservation\tAA_Mouse_Seq\tAA_Mouse_Conservation\t")
    OUTFILE.write("NT_Rhesus_Seq\tNT_Rhesus_Conservation\tAA_Rhesus_Seq\tAA_Rhesus_Conservation\t")
    OUTFILE.write("NT_Marmoset_Seq\tNT_Marmoset_Conservation\tAA_Marmoset_Seq\tAA_Marmoset_Conservation\t")
    OUTFILE.write("NT_GreenMonkey_Seq\tNT_Green Monkey_Conservation\tAA_Green Monkey_Seq\tAA_Green Monkey_Conservation\t")
    OUTFILE.write("NT_Crab-eatingMacaque_Seq\tNT_Crab-eatingMacaque_Conservation\tAA_Crab-eatingMacaque_Seq\tAA_Crab-eatingMacaque_Conservation\t")
    OUTFILE.write("NT_Gibbon_Seq\tNT_Gibbon_Conservation\tAA_Gibbon_Seq\tAA_Gibbon_Conservation\t")
    OUTFILE.write("gnomAD_URL\trsID\tVariantHGVS\tPhenotypeOMIM_info\tPhenotypeList\tCRISPR_Motif_Hits\tNumber_CRISPR_Hits\tNumber_Overlapping_CRISPR_Hits\tClinvarEntrySearch\n")

    # you'll need this hg19.fa file.  I got mine by downloading the chromFA.tar from UCSC, unzipping into the many .fa files, and then concatenating them all
    # This is required for extracting the sequence
    fasta = '%s/hg19.fa'%DataDir

    # This is for the new code that allows monkeys 
    # Updated May 19 2017

    species = { "hg19": "Human",
                "papHam1": "Baboob",
                "otoGar3": "Bushbaby",
                "panTro4": "Chimp",
                "macFas5": "Crab-eating macaque",
                "nomLeu3": "Gibbon",
                "gorGor3": "Gorilla",
                "chlSab1": "Green monkey",
                "calJac3": "Marmoset",
                "ponAbe2": "Orangutan",
                "rheMac3": "Rhesus",
                "saiBol1": "Squirrel monkey",
                "mm10": "Mouse" }
    
    hgsid = conservationWithAA.get_hgsid_from_ucsc_table_browser()
    for line in infile:
        line = line.strip('\n')
        cols = line.split('\t')
        #Extract the variants from the file and put them in a list, this is hard coded, but from column 16-end should only be variants
        # the variants are represented in chr_pos_ref_alt
        variants = cols[15].split('|')
        print variants
        for var in variants:
            if len(var) < 2:
                continue
            print var
            chromosome,position,ref,alt,rsid=var.split('_')
            position = int(position)
            start = position-(WindowSize/2)
            stop = position+(WindowSize/2)
            print "%s\t%d\t%d\n"%(chromosome,start,stop)

#conservation, relies on oriol's conservationWithAA.py program.
            maf = conservationWithAA.get_maf_from_ucsc_table_browser(chromosome,start,stop,hgsid)
            msa = conservationWithAA.get_msa_from_ucsc_table_browser(chromosome,start,stop,hgsid)
            
            # Nucleotide Sequences
            NTsequences,NTidentities=conservationWithAA.get_aligned_nucleotide_sequences_and_percentage_identities_between_human_and_x(maf=maf, x=species.keys())
            AAsequences,AAidentities=conservationWithAA.get_aligned_protein_sequences_and_percentage_identities_between_human_and_x(start,stop,msa=msa, x=species.keys())
#            print species['mm10']
#            print NTsequences['mm10']
#            print NTidentities['mm10']
#            print AAsequences['mm10']
#            print AAidentities['mm10']
#
#            print species['rheMac3']
#            print NTsequences['rheMac3']
#            print NTidentities['rheMac3']
#            print AAsequences['rheMac3']
#            print AAidentities['rheMac3']
#
            #for specie in sorted(NTsequences):
            #    print ("%s\t%s\t%s" % (species[specie], NTsequences[specie], str(NTidentities[specie])))



            #mouse_NT_sequence,mouse_NT_conservation=conservationWithAA.get_nucleotide_percentage_identity_between_human_and_x(maf=maf, x="mm10")
            #try:
            #    AA_conservation=conservationWithAA.get_aminoacid_percentage_identity_between_human_and_x(chromosome, start, stop, maf=maf, msa=msa, x="mm10")
            #except:
            #    AA_conservation=0.
            #    print "Had an error with getting the AA conservation"
            #    print "%s\t%d\t%d"%(chromosome,start,stop)



## gnomAD query.  For now, we'll just copy-paste the gnomAD URL.  Need to automate to get a yes-no for polymorphisms in the region
            gnomAD_url='http://gnomad.broadinstitute.org/region/%s-%d-%d'%(chromosome[3:],start,stop) 

## Fasta sequence from hg19
            seq = pybedtools.BedTool.seq((chromosome,start,stop),fasta)
            #print seq

            # I had to switch to rs ids to capture the indels, since they are represented differently in the variant_summary.txt and the Clinvar file.  Assholes. 
            # So here, I'll check to see if we can use the rs id to link them, if not, then I'll go down the line of evidence to the variant position info, the variant ref and alt alelle, etc.
            if rsid != '.':
                if VarInfoRS.has_key(rsid):
                    rsnum,name,phenos,phenos2 = VarInfoRS[rsid]
                elif VarInfo.has_key('%s_%s_%s_%s'%(chromosome,position,ref,alt)):
                    rsnum,name,phenos,phenos2 = VarInfo['%s_%s_%s_%s'%(chromosome,position,ref,alt)]
                elif VarInfoShort.has_key("%s_%d"%(chromosome,position)):
                    rsnum,name,phenos,phenos2 = VarInfoShort["%s_%d"%(chromosome,position)]
                else:
                    print "no info for this variant"
                    print var
                    rsnum = '.'
                    name = '.'
                    phenos = '.'
                    phenos2 = '.'

            elif VarInfo.has_key('%s_%s_%s_%s'%(chromosome,position,ref,alt)):
                rsnum,name,phenos,phenos2 = VarInfo['%s_%s_%s_%s'%(chromosome,position,ref,alt)]
            elif VarInfoShort.has_key("%s_%d"%(chromosome,position)):
                rsnum,name,phenos,phenos2 = VarInfoShort["%s_%d"%(chromosome,position)]
            else:
                print "no info for this variant"
                print var
                rsnum = '.'
                name = '.'
                phenos = '.'
                phenos2 = '.'

# Motif hits
            PAMsiteHitsFWD,PAMsiteHitsREV = FindMotifInSeq(seq,'NNGRRT')
            #get number of PAM site hits
            PAMsiteHits = PAMsiteHitsFWD + PAMsiteHitsREV
            numhits = len(PAMsiteHits)
            # to look for hits that are within +/- 20bp of the variant, so from positions 130-170 (doesn't matter orientation)
            numCloseHits = 0
            for hit in PAMsiteHits:
                if hit >= 130 and hit <= 170:
                    numCloseHits +=1
            #convert to strings
            PAMsiteHits =map(str,PAMsiteHits)
            hitset = ",".join(PAMsiteHits)

# ClinvarQuery: (12321242[Base Position for Assembly GRCh37]) AND 3[Chromosome] 
            clinvarQuery = '(%d[Base Position for Assembly GRCh37]) AND %s[Chromosome]'%(position,chromosome[3:])
            
# Final write statement
            # Basic stuff
            outfile.write("%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t"%(chromosome,start,stop,position,ref,alt,seq,NTsequences['hg19'],AAsequences['hg19']))
            # conservation stuff
            # Order:  mouse nt seq, mouse nt id, mouse aa seq, mouse aa id
            outfile.write("%s\t%.2f\t%s\t%.2f\t"%(NTsequences['mm10'],NTidentities['mm10'],AAsequences['mm10'],AAidentities['mm10']))
            outfile.write("%s\t%.2f\t%s\t%.2f\t"%(NTsequences['rheMac3'],NTidentities['rheMac3'],AAsequences['rheMac3'],AAidentities['rheMac3']))
            outfile.write("%s\t%.2f\t%s\t%.2f\t"%(NTsequences['calJac3'],NTidentities['calJac3'],AAsequences['calJac3'],AAidentities['calJac3']))
            outfile.write("%s\t%.2f\t%s\t%.2f\t"%(NTsequences['chlSab1'],NTidentities['chlSab1'],AAsequences['chlSab1'],AAidentities['chlSab1']))
            outfile.write("%s\t%.2f\t%s\t%.2f\t"%(NTsequences['macFas5'],NTidentities['macFas5'],AAsequences['macFas5'],AAidentities['macFas5']))
            outfile.write("%s\t%.2f\t%s\t%.2f\t"%(NTsequences['nomLeu3'],NTidentities['nomLeu3'],AAsequences['nomLeu3'],AAidentities['nomLeu3']))
            # clinvar stuff & PAM site hit
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n"%(gnomAD_url,rsnum,name,phenos,phenos2,hitset,numhits,numCloseHits,clinvarQuery))
            
#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    DataDir = '/Users/philliprichmond/Dropbox/Grad School/Gene Therapy/Data'
    infile,outfile,window=GetOptions()
    varinfo,varinfoshort,varinfors= ReadClinvarVariantSummary("%s/Clinvar_VariantSummaryApril2017.txt"%DataDir)
    MakeVarDict(infile,outfile,window,varinfo,varinfoshort,varinfors,DataDir)


