#!/usr/bin/env/python
#This script takes a bam file name, the mean insert size of the bam, and the sample name of the bam, and outputs a config file for use with Pindel
#usage: <script><bam file><insert size><sample name>
import csv
from sys import argv 

script, file1, insertsize, sample = argv
config_name = file1 + "_config.txt"
config = open(config_name, 'w')
config_text = file1 + '\t' + insertsize + '\t' + sample
config.write(config_text)
config.close() 
