import os
import sys
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add
import numpy as np


filepath = "2021-05-26_unmasked.fa"
reference = open("spike_nt.txt", "r").readlines()
reference = list(reference[0].lower())
score = [0]*len(reference)

def mutationchecker(num):
    if reference[num] == sequence[num]:
        return 0
    else:
        if sequence[num] != "n":
            return 1
        else:
            return 0

meta = np.genfromtxt('metadata.csv',delimiter=',',skiprows=1)
print(meta[1,3])

'''
with open(filepath, mode = "r") as handle:
    count = 0
    for record in tqdm(SeqIO.parse(handle, 'fasta')):
        sequence = list(record.seq)
        sequence = sequence[21562:25384]
        try:
            newlist = list(map(mutationchecker, range(0,(len(reference)-1))))
            score = list(map(add, score, newlist))
            count += 1
        except:
            continue

        if count > 10:
            outfile = open("nt_counts.txt", "w")
            for x in score:
                outfile.write(str(x) + "\n")
            outfile.close()
            break
           

print(count)

outfile = open("nt_counts.txt", "w")
for x in score:
    outfile.write(str(x) + "\n")
outfile.close()

print(len(reference))
bases = ["a","c","t","g","n","-"]
otherbases = []
posmutmatrix = np.zeros((len(reference), 6))
with open(filepath, mode = "r") as handle:
    count = 0
    for record in tqdm(SeqIO.parse(handle, 'fasta')):
        sequence = list(record.seq)
        sequence = sequence[21562:25384]
        for nt in range((len(sequence)-1)):
            if sequence[nt] != reference[nt]:
                if sequence[nt] in bases:
                    row = nt
                    col = bases.index(sequence[nt])
                    posmutmatrix[row,col] = posmutmatrix[row,col] + 1
                else:
                    if sequence[nt] in otherbases:
                        continue
                    else:
                        otherbases.append(sequence[nt])
        count += 1
        
        if count > 1000:
            print(otherbases)
            np.savetxt("foo.csv", posmutmatrix, delimiter=",")
            sys.exit()
       

print(otherbases)
np.savetxt("nucleotidemutwise.csv", posmutmatrix, delimiter=",")
'''


