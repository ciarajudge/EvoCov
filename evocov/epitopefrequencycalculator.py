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

def freqcalc(infile, AAepitopes):
    infile = open("Data/"+infile, "r").readlines()
    reference = open("Data/spike_NT.txt", "r").readlines()
    AAreference = open("Data/spike_AA.txt", "r").readlines()
    reference = list(reference[0])
    AAreference = list(AAreference[0])
    NTepitopes = []
    NTwildtypes = []
    epiversions = []
    epiversionscount = []
    bases = []
    posns = []
    for AAepitope in AAepitopes:
        AAepitope = AAepitope.split(",_")[0].split(",")
        AAepitope = [int(x)+1 for x in AAepitope]
        NTepitope = []
        for i in AAepitope:
            NTepitope.append((i*3)-2)
            NTepitope.append((i*3)-1)
            NTepitope.append(i*3)
        NTepitopes.append(NTepitope)
        NTwildtypes.append("".join([reference[x-1] for x in NTepitope]))
        epiversions.append(["".join([reference[x-1] for x in NTepitope])])
        epiversionscount.append([-1])
        bases.append([])
        posns.append([])

    #countries = []
    baseletters = ["A","T","C","G","_"]
    count = 0
    for line in tqdm(infile):
        if list(line)[0] == ">":
            count += 1
            for e in range(0, len(NTwildtypes)):

                
                epi = epitopeadjust(NTwildtypes[e], bases[e], posns[e])
                if epi in epiversions[e]:
                    ind = epiversions[e].index(epi)
                    epiversionscount[e][ind] += 1
                else:
                    epiversions[e].append(epi)
                    epiversionscount[e].append(1)

            #try:
            #    country = line.split("|")[4].split("/")[1].strip()
            #except:
            #    country = line.split("|")[4]
            #if country not in countries:
            #    countries.append(country)
        
            bases = [[] for i in range(len(NTwildtypes))]
            posns = [[] for i in range(len(NTwildtypes))]
        else:
            base = line[-2]
            if base in baseletters:
                position = int(line[1:-2])-1
                for e in range(0, len(NTwildtypes)):
                    if position in NTepitopes[e]:
                        posns[e].append(NTepitopes[e].index(position))
                        bases[e].append(base)

    if not os.isdir("Analysis/Epitopes"):
        os.makedirs("Analysis/Epitopes")
    if not os.isdir("Analysis/Epitopes/CountsandFreqs"):
        os.makedirs("Analysis/Epitopes/CountsandFreqs")
    for e in range(0, len(NTwildtypes)):
        outfile1 = open("Analysis/Epitopes/CountsandFreqs/"+NTwildtypes[e]+"_counts.csv", "w")
        outfile2 = open("Analysis/Epitopes/CountsandFreqs/"+NTwildtypes[e]+"_freqs.csv", "w")
        for i in range(epiversions[e]):
            outfile1.write(epiversions[e][i]+","+str(epiversionscount[e][i])+"\n")
            outfile2.write(epiversions[e][i]+","+str(epiversionscount[e][i]/sum(epiversionscount[e]))+"\n")
            
    countries =  open("Helpers/countries.txt", "r").readlines()
    

freqcalc("difffile.txt", ["448,452,491,492,493,494,495,496,499,_10","352,353,354,355,356,357,358,_3", "502,503,504,505,508,509,510,511,_2"])


sys.exit()


         
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
        if base in bases:
            position = int(line[1:-2])-1
            if position in epitope1:
                posns1.append(epitope1.index(position))
                bases1.append(base)
            elif position in epitope2:
                posns2.append(epitope2.index(position))
                bases2.append(base)
'''

bases1 = []
posns1 = []
bases2 = []
posns2 = []
country = "X"
counting = False
dates = ["2021-06","2021-07"]
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
            date = line.split("|")[1].strip()
            date = "-".join(date.split("-")[0:2])
        except:
            country = line.split("|")[4].split("/")

        if date in dates:
            counting = True
        else:
            counting = False
            
        bases1 = []
        posns1 = []
        bases2 = []
        posns2 = []
    else:
        if counting == False:
            continue
        base = line[-2]
        if base in bases:
            position = int(line[1:-2])-1
            if position in epitope1:
                posns1.append(epitope1.index(position))
                bases1.append(base)
            elif position in epitope2:
                posns2.append(epitope2.index(position))
                bases2.append(base)

np.savetxt("epitope1recent.csv", matrix1, delimiter=',', comments="")
np.savetxt("epitope2recent.csv", matrix2, delimiter=',', comments="")

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

'''




