#!/usr/bin/env python3

'''
Converts a sequence into WYK encoding
based on 2011 paper 'Maximally Efficient Modeling of DNA Sequence Motifs at All
Levels of Complexity'
Usage:
python wyk-encode.py [sequence] <optional flag --diadj>
if --diadj is added it will encode all adjacent dinucleotide, otherwise it will 
just encode mononucleotides
'''

import numpy as np
import argparse

def wyk_encode(base_encoding, seq, adj_dinuc):
    '''
    To Do:
    Given a base_encoding dict and sequence
    If adj_dinuc is false, return list of (base, encoding) for each base
    if adj_dinuc is true, return list of (dinuc, encoding) for each adjacent dinuc
    for dinuc encoding it uses b1 b1xb2 b2 
    '''
    l = []
    for pos, s in enumerate(seq):
        if not adj_dinuc: # mononuc only
            l.append((s, base_encoding[s]))
        else: # encode all adjacent dinucleotides
            di_seq = s + seq[pos+1]
            mono1 = base_encoding[s]
            mono2 = base_encoding[seq[pos+1]]
            outer_prod = np.outer(mono1, mono2).flatten()
            l.append((di_seq, mono1 + list(outer_prod) + mono2))

            if pos == len(seq) - 2: # stop when at last dinucleotide
                break 
    return l

parser = argparse.ArgumentParser()
parser.add_argument('sequence', help='DNA sequence', type=str)
parser.add_argument('-diadj', help='encode adjacent dinucleotides', action='store_true')
args = parser.parse_args()

# wyk encoding:
base_encoding = dict(A=[1,-1,-1], C=[-1,1,-1], G=[-1,-1,1], T=[1,1,1])

encodes = wyk_encode(base_encoding, args.sequence.upper(), args.diadj)

# print out encoding in 2 ways 
# first is the encoding per base/di-base
# second is the encoding for the entire sequence
s = ''
for seq, encode in encodes:
    print(seq, encode)
    e = [str(i) for i in encode]
    s += ','.join(e) + ','

print('\nEntire sequence encoded:')
print(s[:-1]) # remove final comma
