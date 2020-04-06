"""
Pseudocode:
-Def function to turn matrix into a dictionary
-Reading the two sequences
-Def function to compute the score:
    -generating the two matrices F and traceback
    -no need to initialise first row and column, just fill everything with 0s
    -compute the algorithm: SW includes 0 in the max function
    -for each score above thershold, find the score in the whole F matrix and trace back to the 1st 0
-def printing function
-invoke printing function
"""

def mat_to_dict(filepath):
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
    for j in range(len(lists[0])):
        for i in range(len(lists[0])):
            if lists[0][i]+lists[0][j] not in dict:
                dict[lists[0][i]+lists[0][j]]=dict[lists[0][j]+lists[0][i]]
    return dict

s1="HEAGAWGHEE"
s2="PAWHEAE"

def aln_score(s1,s2,d,mat,threshold):
    """This function takes as input 2 strings, the gap penalty, the scoring matrix as a python dict and a threshold
     and returns as a sorted output the alignments with a score above the threshold"""
    M=len(s1)+1
    N=len(s2)+1
    F=[[0]*M for i in range(N)]
    traceback=[[0]*M for i in range(N)]
    for j in range(1,N):
        for i in range(1,M):        #no need to initialise, just apply SW algorithm
            U=F[j-1][i]-d
            L=F[j][i-1]-d
            D=F[j-1][i-1]+mat[s1[i-1]+s2[j-1]]
            F[j][i]=max(U,L,D,0)
            if F[j][i]==U:
                traceback[j][i]="U"
            elif F[j][i]==L:
                traceback[j][i]="L"
            elif F[j][i]==D:
                traceback[j][i]="D"

    best=0
    imax=0
    jmax=0
    results=[]  # This is the list with the complete results
    for j in range(1,N):
        for i in range(1,M):
            if F[j][i]>=threshold:  #it will traceback any alignment with a score>threshold
                best=F[j][i]
                imax=i
                jmax=j
                aln_seq1=""
                aln_seq2=""
                while F[jmax][imax]>0:
                    if traceback[jmax][imax]=="U":
                        aln_seq2+=s2[jmax-1]
                        aln_seq1+="-"
                        jmax-=1
                    elif traceback[jmax][imax]=="L":
                        aln_seq1+=s1[imax-1]
                        aln_seq2+="-"
                        imax-=1
                    elif traceback[jmax][imax]=="D":
                        aln_seq1+=s1[imax-1]
                        aln_seq2+=s2[jmax-1]
                        imax-=1
                        jmax-=1
                par_res=(aln_seq1[::-1],aln_seq2[::-1],"Score: "+str(best))
                results.append(par_res) #this will add each alignment to results

    results.sort(key=lambda tuple: tuple[2], reverse=True) #this sorts results from the highest to the lowest score
    return results

def print_tuple(tuple):
    for i in range(len(tuple)):
        print tuple[i]
    print ""


matrix=mat_to_dict("PAM250_square.txt")
for element in aln_score(s1,s2,2,matrix,20):
    print_tuple(element)