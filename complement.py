# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# AT, CG

def complement(sequence):
    output = ''
    # return the complementary sequence
    for base in sequence:
        if base == 'A':
            output += 'T'
        elif base == 'T':
            output += 'A'
        elif base == 'C':
            output += 'G'
        elif base == 'G':
            output += 'C'
        else:
            output += '*'
    return output
        

testseq = 'ATCGATCGATCGATCN123786'

complement_testseq = complement(testseq)

print(complement_testseq)