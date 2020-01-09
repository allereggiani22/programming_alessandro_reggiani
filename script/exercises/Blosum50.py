from re import *
f= open("./Blosum50.txt","r")
blosum=[]
for line in f:
    line=line.rstrip()
    line=line.split()
    blosum.append(line)
f.close()
for list in blosum[1:]:
    for i in range(1,len(list)):
        if list[i] != "[ARNDCQEGHILKMFPSTWYV]":
            list[i]=int(list[i])

#print blosum

blosum50={}
for i in range(1,len(blosum)):
    for j in range(len(blosum[0])):
        blosum50[blosum[i][0]+blosum[0][j]]=blosum[i][j+1]

#print blosum50["KK"]
        
def score_prot(seq1,seq2):
    if len(seq1)==len(seq2):
        score=0
        for i in range(len(seq1)):
            score = score + blosum50[seq1[i]+seq2[i]]
        return score
    else:
        return "Sequences are different in lenght!"

print score_prot(raw_input("Insert sequence 1: "), raw_input("Insert sequence 2: "))
