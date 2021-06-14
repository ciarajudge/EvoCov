import os
import sys
import subprocess
from tqdm import tqdm
from operator import add
import numpy as np

form = sys.argv[1]
infile = open(sys.argv[2], "r").readlines()
targets = sys.argv[4:]
if form == "NT" :
    reference = open("Data/spike_NT.txt", "r").readlines()
    reference = list(reference[0].lower())
    bases = ["A","C","T","G","N","-"]
    posmutmatrix = np.zeros((len(reference), len(bases), len(targets)))
elif form == "AA" :
    reference = open("Data/spike_AA.txt", "r").readlines()
    reference = list(reference[0].lower())
    bases = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P","S", "T", "W", "Y", "V", "X", "*", "_"]
    posmutmatrix = np.zeros((len(reference)+10, len(bases), len(targets)))
meta = sys.argv[3]
if meta == "location":
    m = 3
elif meta == "variant":
    m = 2
else:
    m = 1

targetcounts = [0]*len(targets)
print(targets)
    
count = 0
otherbases = []
for line in tqdm(infile):
    if list(line)[0] == ">":
        posns = []
        factor = line.split("|")[m].split("\n")[0]
        if any(target == factor for target in targets):
            count += 1
            counting = True
            z = targets.index(factor)
            targetcounts[z] += 1
        else:
            counting = False
    else:
        if counting == True:
            base = line[-2]
            position = int(line[1:-2])
            if position not in posns:
                if base in bases:
                    row = position
                    col = bases.index(base)
                    posmutmatrix[row,col,z] = posmutmatrix[row,col,z] + 1
                    posns.append(position)
                else:
                    if base in otherbases:
                        continue
                    else:
                        otherbases.append(base)
        else:
            continue

print(count)
print(otherbases)
print(targetcounts)
if not os.path.isdir("Analysis/"+meta):
    os.makedirs("Analysis/"+meta)
for x in range(0, len(targets)):
    if targetcounts[x] != 0:
        np.savetxt("Analysis/"+meta+"/"+targets[x]+".csv", posmutmatrix[:,:,x], delimiter=",")
        np.savetxt("Analysis/"+meta+"/"+targets[x]+"frequencies.csv", (posmutmatrix[:,:,x]/targetcounts[x]), delimiter=",")
    else:
        print("No sequences found for "+meta+" "+targets[x])

