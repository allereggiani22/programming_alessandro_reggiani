"""
Pseudocode:
-import s1, s2, BLOSUM62 from input_data.py
-def function that takes as input s1, s2, BLOSUM and gap value to return F (SW) and P matrices
    -defining F and P and filling with 0s, no need to initialise first row/column
    -iterate SW algorithm
    -return F,P
-def function that take F, P, s1 and s2 and returns the best local alignment(s) and corresponding score(s)
    -assign variables to store baest score and relative indexes
    -initialise empty strings to fill during the traceback
    -move backwards from the best score to the first 0 met, using P matrix to insert the correct residue/gap
    -return the two aligned sequences and the score
-def printing function to print results in a nice way

Invoke the 3 functions
"""

#This program is thought for a Python 2.7 interpreter

from input_data import *
a="GAATTC"
b="ACCA"  #I used these for manual calculation, but I saved the program for running with seq1 and seq2 from the module


def score_traceback(s1,s2,matrix,gap):
    """This function takes as input two sequences, a dictionary and a linear gap penalty
    and returns the two matrices of the scores (F) and of traceback (P) in a tuple"""
    M=len(s1)+1     #M and i always refer to s1 (upper side of the matrix)
    N=len(s2)+1     #N and j always refet to s2 (left side of the matrix)
    F=[[0]*M for j in range(N)]
    P=[[0]*M for j in range(N)] #these cmds generate 2 matrices F and P, 0 filled, of N rows(lists) with length M
    for j in range(1,N):
        for i in range(1,M): #from 1 because the first row and colums in a SW matrix are 0-filled
            U=F[j-1][i]-gap
            L=F[j][i-1]-gap
            D=F[j-1][i-1]+matrix[s1[i-1]+s2[j-1]]
            F[j][i]=max(U,L,D,0)    #In these lines I define the 3 possible values and find the max among them and 0,
                                    #since I do not want negative values in the F matrix
            if F[j][i]==D:
                P[j][i]="D"
            elif F[j][i]==U:
                P[j][i]="U"
            elif F[j][i]==L:
                P[j][i]="L"         #Here I fill the traceback matrix with the correct direction to follow
    return F, P         #Remember: this generates a tuple


def local_align(F,P,s1,s2):
    """This function takes as input two matrices of scores and traceback and two sequences,
    finds the best score in the F matrix and traces back from there to the first 0, using the
    traceback matrix. Then returns a tuple with the two aligned sequences and their score"""
    M=len(s1)+1
    N=len(s2)+1
    best=0
    imax=0
    jmax=0
    for j in range(1,N):
        for i in range(1,M):
            if F[j][i]>best:
                best=F[j][i]
                imax=i
                jmax=j      #I want to find the best score in the F matrix and store it and its indexes

    aln_seq1=""
    aln_seq2=""     #Initialisation of the aln sequences that will be filled

    while F[jmax][imax]>0:      #Traceback will end at the first 0 score that is met
        if P[jmax][imax]=="L":
            aln_seq1+=s1[imax-1]
            aln_seq2+="-"       #I add a gap in the sequence towards which I move
            imax-=1             #decrease the used index
        elif P[jmax][imax]=="U":
            aln_seq2+=s2[jmax-1]
            aln_seq1+="-"       #I add a gap in the sequence towards I move
            jmax-=1             #decrease the used index
        elif P[jmax][imax]=="D":
            aln_seq1+=s1[imax-1]
            aln_seq2+=s2[jmax-1]    #no gaps, there is a match
            imax-=1
            jmax-=1                 #decrease both used indexes
    return (aln_seq1[::-1],aln_seq2[::-1],best)


def print_results(tuple):   #this is useful to present results in an organized way
    "This function prints line by line the elements in a tuple with length 3, in which the last is the score"
    for i in range(2):
        print tuple[i]
    print "Score: "+str(tuple[2])

matrices=score_traceback(seq1,seq2,BLOSUM52,2)
aln=local_align(matrices[0],matrices[1],seq1,seq2)
print_results(aln)
