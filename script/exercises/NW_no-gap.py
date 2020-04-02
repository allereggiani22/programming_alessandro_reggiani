"""Pseudocode:
-Definition of function to read the scoring matrix and turn it into a dictionary (possibly triangular)
-Reading the two sequences
-Definition of scoring function(s1,s2,d,mat):
    generation of the two matrices, F and traceback
    initialisation of first row and column to zero (like all the matrix)
    iteration
    finding the best score in the last row/column
    tracing back till the first 0 met in the first line/column
"""

def mat_to_dict(filepath):
    """this functiion takes a triangular matrix as input and returns a dictionary"""
    matrix=open(filepath,"r")
    lists=[]
    for line in matrix:
        line=line.rstrip()
        line=line.split()
        lists.append(line)
    dict={}
    for j in range(1,len(lists)):
        for i in range(len(lists[j])-1):
            dict[lists[j][0]+lists[0][i]]=int(lists[j][i+1])
    for i in range(len(lists[0])):
        for j in range(len(lists[0])):
            if lists[0][i]+lists[0][j] not in dict:
                dict[lists[0][i]+lists[0][j]]=dict[lists[0][j]+lists[0][i]]
    return dict

s1=raw_input("Insert your first DNA sequence here: ")
s2=raw_input("Insert your second DNA sequence here: ")
M=len(s1)+1
N=len(s2)+1


def score_aln(s1,s2,d,mat):
    F=[]
    traceback=[]
    for j in range(N): #generatio and filling in with 0 of the two matrices
        F.append([])
        traceback.append([])
        for i in range(M):
            F[j].append(0)
            traceback[j].append(0)
    for j in range(1,N):
        for i in range(1,M):        #no need to initialise, directly apply the algorithm
            U=F[j-1][i]-d
            L=F[j][i-1]-d
            D=F[j-1][i-1]+mat[s1[i-1]+s2[j-1]]
            F[j][i]=max(U,L,D)
            if F[j][i]==U:
                traceback[j][i]="U"
            elif F[j][i]==L:
                traceback[j][i]="L"
            elif F[j][i]==D:
                traceback[j][i]="D"

    print F
    print traceback

    max_row=max(F[N-1])
    max_col=0
    for j in range(1,N):
        if F[j][-1]>max_col:
            max_col=F[j][-1]
    best=max(max_row,max_col) #finding the highest score in the matrix

    imax=0
    jmax=0
    for j in range(1,N):
        for i in range(1,M):
            if F[j][i]==best:
                jmax=j
                imax=i

    aln_seq1=""
    aln_seq2=""
    while F[jmax][imax]>0:
        if traceback[jmax][imax]=="L":
            aln_seq2+="-"
            aln_seq1+=s1[imax-1]
            imax-=1
        elif traceback[jmax][imax]=="U":
            aln_seq1+="-"
            aln_seq2+=s2[jmax-1]
            jmax-=1
        elif traceback[jmax][imax]=="D":
            aln_seq1+=s1[imax-1]
            aln_seq2+=s2[jmax-1]
            imax-=1
            jmax-=1
    results=(aln_seq1[::-1],aln_seq2[::-1],"Score: "+str(best))
    output=print_tuple(results)
    return output

def print_tuple(tuple):
    for i in range(len(tuple)):
        print tuple[i]





matrix=mat_to_dict("./nucl_mat_triangular.txt")
score_aln(s1,s2,2,matrix)