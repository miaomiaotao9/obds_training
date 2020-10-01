# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 10:34:33 2020

overlap using chromStart chromend

use a bedfile as reference and comparing b bedfile to every chromStart and chromEnd

First attempt: 
    Comparing B (b) start to A (a) chrstart and end
    Comparing B end to A chrstart and end
    
TO DO:
1) Print number of intersects

2) Print overlapping coordinates from one file
    
3) Include arg.parse() etc 

4) zip at the end

5) add log files l...
    
6) Make sure both bed file chr match up 

@author: Jennifer, Hebe, Sophie and Andrés
"""



# a = /Users/jennifer/Desktop/week2/week2/test_overlap/a_dnase1_chr21.bed
# b = /Users/jennifer/Desktop/week2/week2/test_overlap/b_dnase1_chr21.bed

import argparse
import logging as L

parser = argparse.ArgumentParser()
parser.add_argument('--input1', '-i1', dest='abedfilepath', help='input path file')
parser.add_argument('--input2', '-i2', dest='bbedfilepath', help='input path file')
parser.add_argument('--output', '-o', dest='outputfilepath', help='output count path file')
args = parser.parse_args()
# parser.add_argument('--padding', '-p', dest='bedfilepad', type= int, default=0, help='number of bases added')
# args = parser.parse_args()


# with open(args.samfilepath, 'r') as samfile:
    # with gzip.open(args.bedfilepath, 'wt') as bedfile:

overlap_count=0   

with open(args.abedfilepath,'r') as abedfile, open(args.bbedfilepath, 'r') as bbedfile:
    with open(args.outputfilepath, 'wt') as outputfile:
        for line_a in abedfile:
            col_a = line_a.split()
            chromStart_a = int(col_a[1])
            chromEnd_a = int(col_a[2])
            # print(line_a)
            # print(col_a)
            #print(chromStart_a,chromEnd_a)
            #reset iteration of fileb 
            bbedfile.seek(0)
            L.info("start inner loop")
            for line_b in bbedfile:
                col_b = line_b.split()
                chromStart_b = int(col_b[1])
                chromEnd_b = int(col_b[2])
                # print(chromStart_a,chromStart_b, chromEnd_a, chromEnd_b)
                if (chromStart_a >= chromStart_b) and (chromStart_a <= chromEnd_b):
                    overlap_count += 1
                    # print(overlap_count) 
                    outputfile.write(f'{overlap_count}\t{chromStart_a}\t{chromEnd_a}\n')
                elif (chromEnd_a >= chromStart_b) and (chromEnd_a <= chromEnd_b):
                    overlap_count += 1
                    # print(overlap_count)
                    outputfile.write(f'{overlap_count}\t{chromStart_a}\t{chromEnd_a}\n')
                    
                    
#--input1 /Users/Sophie/Documents/UniversityofOxford/DPhil/OBDS/b_dnase1_chr21.bed --input2 /Users/Sophie/Documents/UniversityofOxford/DPhil/OBDS/a_dnase1_chr21.bed -o /Users/Sophie/Documents/UniversityofOxford/DPhil/OBDS/overlapoutput.txt
            
# --input1 '/Users/hchen/desktop/OBDS/Downloads/Course material/week2/2309 Wednesday/test_overlap/brain_dnase1_chr21.bed' --input2 '/Users/hchen/desktop/OBDS/Downloads/Course material/week2/2309 Wednesday/test_overlap/gut_dnase1_chr21.bed' -o '/Users/hchen/desktop/OBDS/Downloads/Course material/week2/2309 Wednesday/test_overlap/overlapoutput.txt'

      
                  
 