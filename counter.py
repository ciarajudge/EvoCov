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
    score = [0]*len(reference)
elif form == "AA" :
    reference = open("Data/spike_AA.txt", "r").readlines()
    reference = list(reference[0].lower())
    score = [0]*(len(reference)+10)
print(len(reference))
count = 0
for line in tqdm(infile):
    if list(line)[0] == ">":
        count += 1
        posns = []
    else:
        position = int(line[1:-2])
        if position not in posns:
            #print(position)
            score[position] += 1
            posns.append(position)

print(count)
outfile = open(sys.argv[3], "w")
for x in score:
    outfile.write(str(x) + "\n")
outfile.close()

