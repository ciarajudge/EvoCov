import sys
import os
import numpy as np
import csv
import math
from tqdm import tqdm
import json
import wget
import requests
import datetime
import subprocess
from itertools import repeat
import operator
from operator import itemgetter

def epitopeadjust(epitope, bases, posns):
    if len(posns) > 0:
        for i in range(0, len(posns)):
            epitope = list(epitope)
            epitope[posns[i]] = bases[i]
            epitope = "".join(epitope)
    return epitope

infile = open("Data/difffile.txt", "r").readlines()
reference = open("Data/spike_NT.txt", "r").readlines()
AAreference = open("Data/spike_AA.txt", "r").readlines()
reference = list(reference[0])
AAreference = list(AAreference[0])
AAepitope1 = [449,453,492,493,494,495,496,497,500]
AAepitope2 = [353,354,355,356,357,358,359]
epitope1 = []
for i in AAepitope1:
    epitope1.append((i*3)-2)
    epitope1.append((i*3)-1)
    epitope1.append(i*3)
epitope2 = []
for i in AAepitope2:
    epitope2.append((i*3)-2)
    epitope2.append((i*3)-1)
    epitope2.append(i*3)    

epitope1seq = [reference[x-1] for x in epitope1]
epitope2seq = [reference[x-1] for x in epitope2]
print("".join(epitope1seq))
print("".join(epitope2seq))

wildtypeepi1 = "".join(epitope1seq)
wildtypeepi2 = "".join(epitope2seq)

epi1versions = [wildtypeepi1]
epi2versions = [wildtypeepi2]
epi1versionscount = [0]
epi2versionscount = [0]
countries = []

bases1 = []
posns1 = []
bases2 = []
posns2 = []
for line in tqdm(infile):
    if list(line)[0] == ">":
        epi1 = epitopeadjust(wildtypeepi1, bases1, posns1)
        epi2 = epitopeadjust(wildtypeepi2, bases2, posns2)
        if epi1 in epi1versions:
            ind = epi1versions.index(epi1)
            epi1versionscount[ind] += 1
        else:
            epi1versions.append(epi1)
            epi1versionscount.append(1)
        if epi2 in epi2versions:
            ind = epi2versions.index(epi2)
            epi2versionscount[ind] += 1
        else:
            epi2versions.append(epi2)
            epi2versionscount.append(1)
        try:
            country = line.split("|")[4].split("/")[1].strip()
        except:
            country = line.split("|")[4]
        if country not in countries:
            countries.append(country)
        
        bases1 = []
        posns1 = []
        bases2 = []
        posns2 = []
    else:
        base = line[-2]
        if base != "N":
            position = int(line[1:-2])-1
            if position in epitope1:
                posns1.append(epitope1.index(position))
                bases1.append(base)
            elif position in epitope2:
                posns2.append(epitope2.index(position))
                bases2.append(base)

print(epi1versions)
print(epi1versionscount)
print(epi2versions)
print(epi2versionscount)
print(countries)

matrix1 = np.zeros((len(epi1versions), len(countries)))
matrix2 = np.zeros((len(epi2versions), len(countries)))
'''
bases1 = []
posns1 = []
bases2 = []
posns2 = []
country = "X"
for line in tqdm(infile):
    if list(line)[0] == ">":
        if country != "X":
            epi1 = epitopeadjust(wildtypeepi1, bases1, posns1)
            epi2 = epitopeadjust(wildtypeepi2, bases2, posns2)
            ind = epi1versions.index(epi1)
            countryind = countries.index(country)
            matrix1[ind, countryind] += 1
            ind = epi2versions.index(epi2)
            matrix2[ind, countryind] += 1

        try:
            country = line.split("|")[4].split("/")[1].strip()
        except:
            country = line.split("|")[4].split("/")
        bases1 = []
        posns1 = []
        bases2 = []
        posns2 = []
    else:
        base = line[-2]
        if base != "N":
            position = int(line[1:-2])-1
            if position in epitope1:
                posns1.append(epitope1.index(position))
                bases1.append(base)
            elif position in epitope2:
                posns2.append(epitope2.index(position))
                bases2.append(base)


np.savetxt("epitope1.csv", matrix1, delimiter=',', comments="")
np.savetxt("epitope2.csv", matrix2, delimiter=',', comments="")
'''

countryoutfile = open("countries.txt", "w")
for x in countries:
    countryoutfile.write(x+"\n")

epitope1outfile = open("epitope1versions.txt", "w")
for x in epi1versions:
    epitope1outfile.write(x+"\n")
    
                
epitope2outfile = open("epitope2versions.txt", "w")
for x in epi2versions:
    epitope2outfile.write(x+"\n")






