import sys, os

# This dirty script takes in 1kGP.3202_samples.pedigree_info.txt 
# which lists the relationships among individuals
# and makes a ped file for each trio.
# this file looks like this:
# sampleID fatherID motherID sex
# HG01711 HG01709 HG01710 2


infile = open('1kGP.3202_samples.pedigree_info.txt','r')
# pop off the header line
headerline = infile.readline()

for line in infile:
    sampleid,fatherid,motherid,sex = line.strip('\n').split(' ')
    if fatherid=='0' or motherid=='0':
        continue
    # now we know we're dealing with a line that has a father+mother, so we build out a new ped file based on the sample id
    newpedfile = open('%s.ped'%sampleid,'w')
    newpedfile.write(headerline)
    newpedfile.write('%s\t0\t0\t1\n'%fatherid)
    newpedfile.write('%s\t0\t0\t2\n'%motherid)
    newpedfile.write('%s\t%s\t%s\t%s\n'%(sampleid,fatherid,motherid,sex))

