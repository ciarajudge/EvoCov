import os
import sys
import time
import subprocess
from tqdm import tqdm
from Bio import Align as align
from operator import add
import numpy

#####Functions#####
def mutationchecker(num):
    if spikeseq[num] == sequence[num]:
        return 0
    else:
        if sequence[num] != "X":
            return 1
        else:
            return 0

def uniquemutationchecker(num):
    if spikeseq[num] == sequence[num]:
        return 0
    else:
        if sequence[num] in mutationdict[num]:
            return 0
        else:
            mutationdict[num].append(sequence[num])
            return 1

def frameshiftfixer(sequence):
    sequence = list(sequence)
    if sequence[68:70] == ["S", "G"]:
        sequence.insert(68, "X")
        sequence.insert(68, "X")
    if sequence[143:145] == ["Y", "H"]:
        sequence.insert(143, "X")
    return sequence

#####Read in Spike Sequence#####
print("Reading in reference Spike Sequence")
spikeseq = open("spike.txt", "r").readlines()
spike = ""
for l in spikeseq:
    spike = spike + l.replace("\n", "")
spikeseq = list(spike)

"""
######Read in Gisaid Sequences#####
print("Parsing Spike Sequences from GISAID")
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

#####Fix frameshift#####
print("Detecting and fixing alignment errors")
sequences = list(map(frameshiftfixer, sequences))

#####Make dictionaries for later#####
countries = [i.split("|")[-1].replace("\n", "") for i in labels]
dates = [i.split("|")[2] for i in labels]
metaandseqdict = [{'date': date, 'country': country, 'sequence': sequence} for date,country,sequence in zip(dates,countries,sequences)]


outfile = open("correctsequences.txt", "w")
for x in sequences:
    outfile.write(str(x) + "\n")
outfile.close()
outfile = open("correctlabels.txt", "w")
for x in labels:
    outfile.write(str(x) + "\n")
outfile.close()
#numpy.save('sequencedict.npy', metaandseqdict) 
"""
#####Read in Sequence Dict#####
"""
print("Reading in Sequence Dictionary")
metaandseqdict = numpy.load('sequencedict.npy',allow_pickle='TRUE').item()
print("Sequence Dictionary Loaded")
sequences = metaandseqdict["sequence"]
len(sequences)
"""
#####Simple Counts#####
"""
count = 0
for sequence in tqdm(sequences):
    try:
        newlist = list(map(mutationchecker, range(0,(len(spikeseq)-1))))
        score = list(map(add, score, newlist))
        count += 1
    except:
        continue


outfile = open("fixedcounts.txt", "w")
for x in score:
    outfile.write(str(x) + "\n")
outfile.close()

print(count)
sys.exit()
"""


#####Unique Mutations Count#####
"""
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
outfile = open("fixeduniquecounts.txt", "w")
for x in score:
    outfile.write(str(x) + "\n")
outfile.close()
"""


#####Mutation Tally#####
"""
muttable = numpy.zeros(21, 222)
for sequence in tqdm(sequences):
    try:
        newlist = list(map(uniquemutationchecker, range(0,(len(spikeseq)-1))))
        score = list(map(add, score, newlist))
        count += 1
    except:
        continue

"""


#####Temporal Analysis#####
sequences = open("correctsequences.txt", "r").readlines()

infile = open("correctlabels.txt", "r")
labels = []
for line in infile:
        if line.strip():
                labels.append(line)
        else:
                continue

dates = [i.split("|")[2] for i in labels]
years = [i.split("-")[0] for i in dates]


_2019score = [0]*len(spikeseq)
_2020score = [0]*len(spikeseq)
_2021score = [0]*len(spikeseq)
_2019count = 0
_2020count = 0
_2021count = 0
count = -1
for sequence in tqdm(sequences):
    count += 1
    sequence = sequence.replace("['", "").replace("', '","").replace("']","")
    try:
        newlist = list(map(mutationchecker, range(0,(len(spikeseq)-1))))
        if years[count] == '2019':
            _2019score = list(map(add, _2019score, newlist))
            _2019count += 1
        elif years[count] == '2020':
            _2020score = list(map(add, _2020score, newlist))
            _2020count += 1
        else:
            _2021score = list(map(add, _2021score, newlist))
            _2021count += 1
    except:
        continue


outfile = open("2019counts.txt", "w")
for x in _2019score:
    outfile.write(str(x) + "\n")
outfile.close()

outfile = open("2020counts.txt", "w")
for x in _2020score:
    outfile.write(str(x) + "\n")
outfile.close()

outfile = open("2021counts.txt", "w")
for x in _2021score:
    outfile.write(str(x) + "\n")
outfile.close()

print(_2019count)
print(_2020count)
print(_2021count)
sys.exit()




