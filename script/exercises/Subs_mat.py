#Compute the substitution score matrix starting from these alignments
from math import *

seq1="ACAGGTGGACCTCTATATGG"
seq2="ACTGGTCGACTTCCGGATCG"

def freq_nucl(seq1,se12):
    freq_sing={}
    nucl=["A","C","T","G"]
    for base in nucl:
        freq_sing[base]=((seq1+seq2).count(base))/float(len(seq1+seq2))
    return freq_sing

def freq_comb(seq1,seq2):
    if len(seq1)==len(seq2):
        nucl=["A","C","T","G"]
        freq_pair={}
        fn=freq_nucl(seq1,seq2)
        for i in nucl:
            for j in nucl:
                freq_pair[i+j]=1
        for n in range(len(seq1)):
            freq_pair[seq1[n]+seq2[n]]=freq_pair.get(seq1[n]+seq2[n],1)+1
        for key in freq_pair:
            freq_pair[key]=log((freq_pair[key]/float(len(seq1)))/(fn[key[0]]*fn[key[1]]))
        return freq_pair
    else:
        return "Sequences are different in lenght!"

print freq_comb(seq1,seq2)

