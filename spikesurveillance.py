import os
import sys
import time
import subprocess
from tqdm import tqdm
from Bio import Align as align
from operator import add
import numpy

spikeseq = open("spike.txt", "r").readlines()
spike = ""
for l in spikeseq:
    spike = spike + l.replace("\n", "")
spikeseq = list(spike)

sequencefile = open(sys.argv[1], "r").readlines()
labels = []
sequences = []
score = [0]*len(spikeseq)
s = ""
for l in sequencefile:
    if ">" in l:
        labels.append(l)
        sequences.append(s.replace("\n", ""))
        s = ""
    else:
        s = s+l
sequences = sequences[1:len(sequences)]

def mutationchecker(num):
    if spikeseq[num] == sequence[num]:
        return 0
    else:
        if sequence[num] != "X":
            return 1
        else:
            return 0

countries = [i.split("|")[-1].replace("\n", "") for i in labels]
dates = [i.split("|")[2] for i in labels]

metaandseqdict = [{'date': date, 'country': country, 'sequence': sequence} for date,country,sequence in zip(dates,countries,sequences)]

print(metaandseqdict[1])
"""
count = 0
for sequence in tqdm(sequences):
    try:
        newlist = list(map(mutationchecker, range(0,(len(spikeseq)-1))))
        score = list(map(add, score, newlist))
        count += 1
    except:
        continue


outfile = open("counts.txt", "w")
for x in score:
    outfile.write(str(x) + "\n")
outfile.close()
"""

def uniquemutationchecker(num):
    if spikeseq[num] == sequence[num]:
        return 0
    else:
        if sequence[num] in mutationdict[num]:
            return 0
        else:
            mutationdict[num].append(sequence[num])
            return 1



keys = range(0, 1273)
mutationdict = {key: ["X"] for key in keys}

count = 0
for sequence in tqdm(sequences):
    try:
        newlist = list(map(uniquemutationchecker, range(0,(len(spikeseq)-1))))
        score = list(map(add, score, newlist))
        count += 1
    except:
        continue

print(count)
outfile = open("uniquecounts.txt", "w")
for x in score:
    outfile.write(str(x) + "\n")
outfile.close()
