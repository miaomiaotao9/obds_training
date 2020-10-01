# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:18:45 2020

Bamtobed pysam
# BED file format 
1) chromosome 
2) chrom start (0 based)) 
3) chrom end 
4) Name 
5) Score 
6) Strand 

@author: whole OBDS team

To do:
    1) v/ Finish BAM to BED file conversion script
    2) v/ Make command line runnable
    3) Make it read from stdin and write to stdout
    4) Add an option to write to a compressed file
    5) Add an option to truncate read coordinates to the first base
"""
# Convert Sam to BED using Pysam
# read in Sam file in Pysam
# iterate over the sam file line by line
# Create the bed file
# import pysam

import pysam
import argparse
import gzip
import sys 

parser = argparse.ArgumentParser(prog="This script does BAM to BED conversion")
parser.add_argument('--input', '-i', dest='bamfilepath', help='Input bam file path')
parser.add_argument('--output', '-o', dest='outputfilepath', help='Output bed file path')
parser.add_argument('--gzip','-z', dest = 'gzipbedfile', action= "store_true", default=False, help = 'Compressed the output bed files')
args = parser.parse_args()

# outputfilepath = "/Users/andresnoe/obds_sep20/working_directory/outputfiletest.bed"
# bamfilepath = "/Users/andresnoe/obds_sep20/working_directory/ERR1755082.test.sort.bam"



bamfile = pysam.AlignmentFile(args.bamfilepath, "rb")

if args.gzipbedfile:
    outputfile = gzip.open(args.outputfilepath, 'wt')
elif args.outputfilepath == "-":
    outputfile = sys.stdout
else:
    outputfile = open(args.outputfilepath, 'w')
    #when put bam file in stdin - do not have index 
for alignment in bamfile.fetch(until_eof=True):
            # pysam iterates over header differently and puts it in a 'different place'
    if alignment.is_paired:
        outputfile.write(f'''{alignment.reference_name}\t{alignment.reference_start}\t{alignment.reference_end}\t{alignment.query_name}\t{alignment.mapping_quality}\t.\n''')

outputfile.close()
bamfile.close()
        

