#! /usr/bin/env python

"""
Function:
    Compare BED files 1 and 2,
    For any intervals in 1 that overlaps with 2
    Output that interval_1 (ie comparing 1 to 2 is not the same as comparing 2 to 1)
    
Input:
    BED files 1 and 2
Output:
    BED files of intervals_1 that overlaps with 2

Authors:
    The A Team
    Harry
    Hongyan
    Maisha
    Phil
"""

# read BED 1, for each line in BED 1, read BED 2 to find overlapping region(s)
# during development, we used gut as BED 1, and brain as BED 2


# input files
# bedfilepath_1= '/Users/mjabeen/Documents/DPhil_2020/OBDS/obds_wd/overlapbed/gut_dnase1_chr21.bed'
# bedfilepath_2 = '/Users/mjabeen/Documents/DPhil_2020/OBDS/obds_wd/overlapbed/brain_dnase1_chr21.bed'
# output_overlapbed = '/Users/mjabeen/Documents/DPhil_2020/OBDS/obds_wd/overlapbed/dnase1_chr21_overlap.bed'


import argparse 
import logging as L


parser = argparse.ArgumentParser()
parser.add_argument('--input1', '-1', dest='bedfile1', help='input1 file path')
parser.add_argument('--input2', '-2', dest='bedfile2', help='input2 file path')
parser.add_argument('--output','-o', dest='overlap', help='output')
parser.add_argument('--log', dest='logconfig', default='INFO', help='verbosity of logging')
args=parser.parse_args()

L.basicConfig(level=args.logconfig)
L.info(f"Parsed arguments as {repr(args)}")

# data structure of BED files we work with:
#     col 1 = chromosome
#     col 2 = start
#     col 3 = end
#     col 4 = score

def getstart(line):
    fields = line.split()
    startpos = int(fields[1])
    return startpos

def getend(line):
    fields = line.split()
    endpos = int(fields[2])
    return endpos
 
def getchrom(line):
    fields = line.split()
    chrom = fields[0]
    return chrom

L.info("Initializing the counter for overlapping regions.")    
overlap_count = 0  
line_1_count = 0


L.info(f"Opening BED file 1 at: {args.bedfile1}")
with open(args.bedfile1, 'r') as bedfile_1:
    L.info(f"Creating output file at: {args.overlap}")
    with open (args.overlap, 'wt') as output:
        for line_1 in bedfile_1: # read line by line
            line_1_count += 1
            startpos_1 = getstart(line_1)
            endpos_1 = getend(line_1)
            chrom = getchrom(line_1)    
            line_comparison = 0
            
            L.debug(f"Working on line {line_1_count} of BED file 1")
            L.debug("Opening BED file 2")
            with open(args.bedfile2, 'r') as bedfile_2:    
                line_2_count = 0
                for line_2 in bedfile_2:
                    line_2_count += 1
                    startpos_2 = getstart(line_2)
                    endpos_2 = getend(line_2)
                    
                    L.debug(f"Checking if BED file 2 overlaps with BED file 1 at {line_2_count}; iteration number {line_1_count}")
                    if startpos_2 in range(startpos_1, endpos_1 + 1) or endpos_2 in range(startpos_1, endpos_1 + 1):
                        overlap_count += 1
                        line_comparison += 1   
        
            L.debug("Writing into output file")
            output.write(f'{chrom}\t{startpos_1}\t{endpos_1}\t{line_comparison}\n')
        
    L.info(f'The total number of overlapping intervals of BED file 1 in BED file 2 is: {overlap_count}')
                
""" Output file format:
    - no. overlap per line in bedfile1
    - which line inserted output file
    - Bedfile:
            chrm    
            chrm start pos
            chrm end
            overlapping intervals """
            
            
    
                


























