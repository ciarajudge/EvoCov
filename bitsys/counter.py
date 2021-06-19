import os
import sys
import subprocess
from tqdm import tqdm
from operator import add
import numpy as np

form = sys.argv[1]
infile = open(sys.argv[2], "r").readlines()
if form == "NT" :
    reference = open("Data/spike_NT.txt", "r").readlines()
    reference = list(reference[0].lower())
    bases = ["A","C","T","G","N","-"]
    posmutmatrix = np.zeros((len(reference), len(bases)))
elif form == "AA" :
    reference = open("Data/spike_AA.txt", "r").readlines()
    reference = list(reference[0].lower())
    bases = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P","S", "T", "W", "Y", "V", "X", "*", "_"]
    posmutmatrix = np.zeros((len(reference)+10, len(bases)))
outfile = sys.argv[3]

print(len(reference))

count = 0
otherbases = []
for line in tqdm(infile):
    if list(line)[0] == ">":
        count += 1
        posns = []
    else:
        base = line[-2]
        position = int(line[1:-2])
        if position not in posns:
            if base in bases:
                row = position
                col = bases.index(base)
                posmutmatrix[row,col] = posmutmatrix[row,col] + 1
                posns.append(position)
            else:
                if base in otherbases:
                    continue
                else:
                    otherbases.append(base)

print(count)
print(otherbases)
np.savetxt(outfile, posmutmatrix, delimiter=",")