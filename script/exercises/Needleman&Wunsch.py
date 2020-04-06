"""pseudocode:
- function to convert substitution matrix into dictionary
- read (from file) sequences 1 and 2
- scoring function definition(seq1,seq2,d,matrix):
    - generation of 2 matrices (F and traceback)
    - initialisation of first row and column
    - application of the algorithm to iterate
    - reconstruction of the traceback
    - returning aln_seq1, aln_seq2 and score
- definition of a printing function to present results properly
- calling the 3 functions
"""

def mat_to_dict(filepath):
    """This function returns a dictionary from a space separated square matrix, given the filepath"""
    matrix=open(filepath, "r")
    dict={}
    lines=[]
    for line in matrix:
        line=line.rstrip()
        line=line.split()
        lines.append(line)
    for i in range(len(lines[0])):
        for j in range(1,len(lines)):
            dict[lines[0][i]+lines[j][0]]=int(lines[j][i+1])
    matrix.close
    return dict

s1=raw_input("Enter your first DNA sequence here: ")
s2=raw_input("Enter your second DNA sequence here: ")

def global_align(s1,s2,d,mat):
    F=[]
    traceback=[]
    M=len(s1)+1
    N=len(s2)+1
    for j in range(N): # j = ROW = s2
        F.append([])
        traceback.append([])
        for i in range(M):  # i = COLUMN = s1
            F[j].append(0)
            traceback[j].append(0)
    for i in range(1,len(F[0])):        #initialise first row
        F[0][i]=F[0][i-1]-d
        traceback[0][i] = "L"
    for j in range(1,len(F)):       #initialise first column
        F[j][0]=F[j-1][0]-d
        traceback[j][0]="U"
    for i in range(1,M):
        for j in range(1,N):   #iteration
            U=F[j-1][i]-d
            L=F[j][i-1]-d
            D=F[j-1][i-1]+mat[s1[i-1]+s2[j-1]]
            F[j][i]=max(U,L,D)
            if F[j][i]==U:              #filling traceback matrix
                traceback[j][i]="U"
            elif F[j][i]==L:
                traceback[j][i]="L"
            elif F[j][i]==D:
                traceback[j][i]="D"     #ok till here


  #now I must rebuild the two aligned sequences, following the traceback, until traceback[0][0]

    imax=M-1
    jmax=N-1
    aln_seq1=""
    aln_seq2=""

    while imax>0 and jmax>0:                 #if I move up, I put the gap in seq1
        if traceback[jmax][imax]=="U":
            aln_seq1+="-"
            aln_seq2+=s2[jmax-1]
            jmax-=1
        elif traceback[jmax][imax]=="L":    #if I move left, I put the gap in seq2
            aln_seq2+="-"
            aln_seq1+=s1[imax-1]
            imax-=1
        elif traceback[jmax][imax]=="D":    #if I move on the diagonal, I allow the match/mismatch
            aln_seq1+=s1[imax-1]
            aln_seq2+=s2[jmax-1]
            imax-=1
            jmax-=1

    seq1=aln_seq1[::-1]
    seq2=aln_seq2[::-1]
    max_score="Score: "+ str(F[N-1][M-1])

    results=(seq1,seq2,max_score)
    return results

def print_tuple(tuple):
    """This fuction prints line by line the content of a tuple"""
    for i in range(len(tuple)):
        print tuple[i]

mat=mat_to_dict("nucl_matrix.txt")
print_tuple(global_align(s1,s2,2,mat))
