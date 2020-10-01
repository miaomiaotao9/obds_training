#3
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