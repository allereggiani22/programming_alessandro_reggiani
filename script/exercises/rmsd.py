file1=open("./protein1.pdb","r")
file2=open("./protein2.pdb","r")
from math import *
model1=[]
model2=[]
for line in file1:
	line=line.rstrip()
	if "CA" in line:
		list1=line.split()
		model1.append([float(list1[6]),float(list1[7]),float(list1[8])])
for line in file2:
	line=line.rstrip()
	if "CA" in line:
		list2=line.split()
		model2.append([float(list2[5]),float(list2[6]),float(list2[7])])
#this first part generates two lists of lists, one for each file,
#in which every element represents the coordinates of each alfa carbon
def rmsd(a,b):
	"""This function computes the RMSD of two sequences, given the two lists of alfa carbons coordinates"""
	Di=[]
	for i in range(len(a)):
		D=((a[i][0]-b[i][0])**2+(a[i][1]-b[i][1])**2+(a[i][2]-b[i][2])**2)
		Di.append(D)
	RMSD=sqrt(0.01*fsum(Di))
	return RMSD

print rmsd(model1,model2)
file1.close()
file2.close()
