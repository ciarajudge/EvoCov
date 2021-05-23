import os
import sys
import time
import subprocess
from tqdm import tqdm
from Bio import Align as align
from operator import add

spikeseq = open("spike.txt", "r").readlines()
spike = ""
for l in spikeseq:
    spike = spike + l.replace("\n", "")
spikeseq = list(spike)
print(spikeseq)

sequencefile = open(sys.argv[1], "r").readlines()
labels = []
sequences = []
score = []*len(spikeseq)
s = ""
for l in sequencefile:
    if ">" in l:
        labels.append(l)
        sequences.append(s.replace("\n", ""))
        s = ""
    else:
        s = s+l
sequences = sequences[1:]
print(len(sequences))

def mutationchecker(num):
    if spikeseq[num] == sequence[num]:
        return 0
    else:
        if sequence[num] != "N":
            return 1
        else:
            return 0

for sequence in tqdm(sequences):
    try:
        newlist = list(map(mutationchecker, range(0,(len(spikeseq)-1))))
        score = list(map(add, score, newlist))
    except:
        continue

print(score)
outfile = open("counts.txt", "w")
for x in score:
    outfile.write(score + "\n")
outfile.close()

