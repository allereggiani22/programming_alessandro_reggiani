"""
Define a function to generate dictionary from matrix
do the same for other matrix
read sequences from fasta file
define function with arguments (seq1, seq2, matrix) to score the alignments
call the function for each matrix and each alignment
print results in some nice way

"""

def matrix_to_dict(filepath):
    file=open(filepath,"r")
    matrix={}
    aacids=[]
    for line in file:
        line=line.rstrip()
        line=line.split()
        aacids.append(line)
    for i in range(len(aacids[0])):
        for j in range(1,len(aacids)):
            matrix[aacids[0][i]+aacids[j][0]]=int(aacids[i+1][j])
    file.close()
    return matrix

BLOSUM=matrix_to_dict("./BLOSUM62_square.txt")
PAM=matrix_to_dict("./PAM250_square.txt")

aln1=[]
fasta=open("./alignments.fasta","r")
for line in fasta:
    if not line.startswith(">"):
        aln1.append(line.rstrip())
fasta.close()

aln2=[aln1[2],aln1[3]]
aln3=[aln1[4],aln1[5]]
aln1=aln1[0:2]

#print aln1
#print aln2
#print aln3

aln=[aln1,aln2,aln3]

def score(seq1,seq2,matrix):
    if len(seq1) == len(seq2):
        aln_seq=seq1+"\n"+seq2
        score=0
        for i in range(len(seq1)):
            score=score+matrix.get(seq1[i]+seq2[i],-2)
    return [aln_seq,score,]

def print_result(list):
    print list[0],"\nScore:", list[1]


PAM_scores=[]
BLOSUM_scores=[]
for list in aln:
    p=score(list[0],list[1],PAM)
    PAM_scores.append(p)
    b=score(list[0],list[1],BLOSUM)
    BLOSUM_scores.append(b)

print "PAM scores:\n"
for element in PAM_scores:
    print_result(element)
print " "
print "BLOSUM scores:\n"
for element in BLOSUM_scores:
    print_result(element)
