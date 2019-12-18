seq1=raw_input("insert seq1: ")
seq2=raw_input("insert seq2: ")
def al_score(seq1,seq2):
    score=0
    for i in range(len(seq1)):
        if seq1[i]==seq2[i]:
            score+=1
        else:
            score-=0
    return score
if len(seq1)!=len(seq2):
    print "Sequences are different in lenght!"
else:
    print al_score(seq1,seq2)