#pseudocode
"""
import math
function to generate a dictionary from matrix file (gaps included)
function to compute a score for each possible alignment, keeping only the best one
    Possible alignments are len(seq1)+len(seq2), and you move one seq on the other by moving the last "-" of one seq
    at the beginning (until seq ends with "-") or the first of the other to the end
    score must be stored ad replaced if new score is higher
    function must return a list with best_seq1 \n best_seq2 as l[0] and best_score as l[1]
function that prints l[0], "\n Score:", l[1]
"""

from math import *

def matrix_to_dict(filepath):
    file=open((filepath),"r")
    nucl_mat={}
    nucl=[]
    for line in file:
        line=line.rstrip()
        line=line.split()
        nucl.append(line)
    for i in range(len(nucl[0])):
        for j in range(1,len(nucl)):
            nucl_mat[nucl[0][i]+nucl[j][0]]=int(nucl[i+1][j])
    file.close()
    return nucl_mat

#print matrix_to_dict("./ungap-matrix.txt")

def ungap_align(seq1,seq2):
    aln_lenght=len(seq1)+len(seq2)
    gaps1=list(len(seq2)*"-")
    gaps2=list(len(seq1)*"-")
    aln=[gaps1+list(seq1),list(seq2)+gaps2]
    max_score=(-2)*len(seq1)
    matrix=matrix_to_dict("./ungap-matrix.txt")
    for i in range(aln_lenght):
        aln_seq1="".join(aln[0])
        aln_seq2="".join(aln[1])
        aln_score=0
        aln_score=aln_score+matrix[aln_seq1[i]+aln_seq2[i]]
        if aln_score>max_score:
            max_score=aln_score
            best_aln=aln_seq1+"\n"+aln_seq2
        if aln[1][-1]=="-":
            aln[1]=aln[1][-1:]+aln[1][:-1]
        else:
            aln[0]=aln[0][1:]+aln[0][:1]
    return [best_aln, max_score]

def print_best(list):
    print list[0],"\n Score:", list[1]

al=ungap_align("AT","CTA")
print_best(al)
#print matrix_to_dict(open("./ungap-matrix.txt","r"))
