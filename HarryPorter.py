# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""






def reverse_complement (sequence):
    output = ''
    for base in sequence:
        if base == 'A':
            output += 'T'
        elif base == 'T':
            output += 'A'
        elif base == 'G':
            output += 'C'
        elif base == 'C':
            output += 'G'
        else:
            print("Unknown base")
    return output[::-1]

test = 'ATGCCCGTTAC123'
reverse_complement(test) 
print(reverse_complement(test))

#

participant_list = ['Andres', 'Hongyan', 'Harry', 'Maisha','Sophie','Phil','Jen']

participant_list.append("David")
participant_list.append("Charlie")
participant_list.append("Kevin")
participant_list.append("Lucy")

print (participant_list)

print (participant_list[2])
print (participant_list[4])

participant_list.sort()
print(participant_list)
print(participant_list[2:5])

participant_dictionary = {}
for participant in participant_list:
    participant_dictionary[participant]= "participan"
 
trainer_list = ["Kevin", "Charlie", "David", "Lucy"]   
 
print(participant_dictionary)
    
for participant in participant_dictionary.values():
    print(participant)

for participant in participant_dictionary.items():
    print(participant)    
    
for participant in participant_dictionary.keys():
    if participant in trainer_list:
        print(participant)
    else:
        continue
    
print(participant_list[0:2])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    