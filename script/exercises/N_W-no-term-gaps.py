"""
Pseudocode:
- Read (from file) s1 and s2
- Read substitution matrix
- Build a "zero filled" matrix ((n+1)*(m+1))
-   intialise F(0,0)
    initialise F[i,0]
    initialise F[0,j]
    iteration --> apply algorithm
    return F
- Fill the first column with gap penalties -(i*d)
- Fill the first row with gap penalties -(j*d)
- For the other positions, use F(i,j) where F(i,j) is maximum[val1,val2,val3)
- Print F in a quite aesthetical way

"""
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


s1="ACTGG"
s2="ACCA"
d=2
subs_mat=matrix_to_dict("nucl_matrix.txt")

def score_m(s1,s2,d,matrix):
    """ Calculates the global dynamic programming (scoring) matrix"""
    F=[] #nice way to generate lists: list comprehension F[[0]*N for x in range(M)]
    traceback=[]
    M=len(s1)+1
    N=len(s2)+1
    for i in range(N):
        F.append([])
        traceback.append([])
        for j in range(M):
            F[i].append(0)
            traceback[i].append(0)
    for i in range(1,N):
        for j in range(1,M):            #iteration
            U=F[i-1][j]-d
            L=F[i][j-1]-d
            D=F[i-1][j-1]+matrix[s1[j-1]+s2[i-1]]
            F[i][j]=max(L,U,D)
            if F[i][j]==D:
                traceback[i][j]="D"
            elif F[i][j]==L:
                traceback[i][j]="L"
            elif F[i][j]==U:
                traceback[i][j]="U"


    max_row=max(F[N-1])
    column=[]
    for i in range(N):
        column.append(F[i][M-1])
    max_col=max(column)
    best=max(max_row,max_col)
    if best in F[N-1]:
        i_max=N-1
        j_max=F[N-1].index(best)
    else:
        i_max=M-1
        j_max=column.index(best)

    # Generates two blank strings to be filled while the aln is built
    # if I move to the left, I keep the nucl from seq1 and add a gap in seq2
    # if I move up, I keep nucl from seq2 and add a gap in seq1
    # alns will be printed reversed at the end

    aln1=""
    aln2=""
    j=j_max
    i=i_max
    while F[i][j]>0:
        if traceback[i][j]=="D":
            aln1+=s1[j-1]
            aln2+=s2[i-1]
            i-=1
            j-=1
        elif traceback[i][j]=="L":
            aln1+=s1[j-1]
            aln2+="-"
            j-=1
        elif traceback[j][i]=="U":
            aln1+="-"
            aln2+=s2[i - 1]
            i-=1
    return aln1[::-1],aln2[::-1],best


def print_aln(tuple):
    """This function prints the result of the scoring function,
    which is a tuple with aln1, aln2 and the score"""
    print tuple[0]
    print tuple[1]
    print "Score:", tuple[2]


alignment=score_m(s1,s2,2,subs_mat)
print_aln(alignment)