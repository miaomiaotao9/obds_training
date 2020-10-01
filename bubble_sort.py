#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 10:29:54 2020

@author: mjabeen
"""

"""Find the maximum number
number_list = [26, 54, 93, 17, 77, 31, 44, 55, 20]
max_number = number_list[0]

Create a loop to find this 
for number in number_list:
    if number > max_number:
        max_number = number
        
print (max_number)"""


"""Sorting - using the selection sort algorithm >>> selection_sort.py"""

"""start by turning above code into function"""

number_list = [26, 54, 93, 17, 77, 31, 44, 55, 20]


# Outer loop - use j to define length of list, going down to zero, excluding last variable
#'Stepping loop'
for j in range (len(number_list), 1, -1):
    for i in range (0, j-1, 1): #Inner loop - 'swapping loop': pushing largest number to the end
        if number_list[i]>number_list[i+1]: #pairwise comparisons
            temp = number_list[i] #creating temporary variable
            number_list [i] = number_list[i+1] #swapping variables (line 31 - 32) 
            number_list[i+1] = temp 
        print(number_list)
        
#instead of i or j think about using meaningful variable names for ease of interpreting code
# e.g. j = step/iteration, i = current index of position in list
    
    
    
    
    
    
    
    
    
    
    
    
    
    