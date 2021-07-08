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


def freqcalc(infile, AAepitopes, dates):
    def epitopeadjust(epitope, bases, posns):
        if len(posns) > 0:
            for i in range(0, len(posns)):
                epitope = list(epitope)
                epitope[posns[i]] = bases[i]
                epitope = "".join(epitope)
        return epitope
    infile = open("Data/"+infile, "r").readlines()
    reference = open("Data/spike_NT.txt", "r").readlines()
    AAreference = open("Data/spike_AA.txt", "r").readlines()
    reference = list(reference[0])
    AAreference = list(AAreference[0])
    NTepitopes = []
    NTwildtypes = []
    AAwildtypes = []
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
        AAwildtypes.append("".join([AAreference[x-1] for x in AAepitope]))
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
    for e in range(0, len(NTwildtypes)): 
        epi = epitopeadjust(NTwildtypes[e], bases[e], posns[e])
        if epi in epiversions[e]:
            ind = epiversions[e].index(epi)
            epiversionscount[e][ind] += 1
        else:
            epiversions[e].append(epi)
            epiversionscount[e].append(1)
    bases = [[] for i in range(len(NTwildtypes))]
    posns = [[] for i in range(len(NTwildtypes))]
    
    if not os.path.isdir("Analysis/Epitopes"):
        os.makedirs("Analysis/Epitopes")
    if not os.path.isdir("Analysis/Epitopes/CountsandFreqs"):
        os.makedirs("Analysis/Epitopes/CountsandFreqs")
    for e in range(0, len(NTwildtypes)):
        outfile1 = open("Analysis/Epitopes/CountsandFreqs/"+AAwildtypes[e]+"_counts.csv", "w")
        outfile2 = open("Analysis/Epitopes/CountsandFreqs/"+AAwildtypes[e]+"_freqs.csv", "w")
        for i in range(len(epiversions[e])):
            outfile1.write(epiversions[e][i]+","+str(epiversionscount[e][i])+"\n")
            outfile2.write(epiversions[e][i]+","+str(epiversionscount[e][i]/sum(epiversionscount[e]))+"\n")
            
    countries =  open("Helpers/countries.txt", "r").readlines()
    countries = [x.strip("\n") for x in countries]
    countrind = 0
    epiversionscount = [[[0 for c in range(len(countries))] for j in range(len(epiversions[i]))] for i in range(len(NTwildtypes))]
    counting = False
    for line in tqdm(infile):
        if list(line)[0] == ">":
            count += 1
            if counting == True:
                for e in range(0, len(NTwildtypes)):
                    epi = epitopeadjust(NTwildtypes[e], bases[e], posns[e])
                    ind = epiversions[e].index(epi)
                    epiversionscount[e][ind][countrind] += 1

            date = line.split("|")[1].strip()
            date = "-".join(date.split("-")[0:2])
            if date in dates:
                counting = True
            else:
                counting = False
            
            try:
                country = line.split("|")[4].split("/")[1].strip()
            except:
                country = line.split("|")[4].strip("\n")

            if country in countries:
                countrind = countries.index(country)
            else:
                counting = False
                '''
                countries.append(country)
                for e in range(len(NTwildtypes)):
                    for c in range(len(epiversions[e])):
                        epiversionscount[e][c].append([])
                countrind = countries.index(country)
                '''
            bases = [[] for i in range(len(NTwildtypes))]
            posns = [[] for i in range(len(NTwildtypes))]
        else:
            if counting == False:
                continue
            base = line[-2]
            if base in baseletters:
                position = int(line[1:-2])-1
                for e in range(0, len(NTwildtypes)):
                    if position in NTepitopes[e]:
                        posns[e].append(NTepitopes[e].index(position))
                        bases[e].append(base)

    if not os.path.isdir("Analysis/Epitopes/ByCountry"):
        os.makedirs("Analysis/Epitopes/ByCountry")
    for e in range(0, len(NTwildtypes)):
        outfile1 = open("Analysis/Epitopes/ByCountry/"+AAwildtypes[e]+".csv", "w")
        for i in range(len(epiversions[e])):
            outfile1.write(epiversions[e][i]+",")
            for c in range(len(countries)):
                outfile1.write(str(epiversionscount[e][i][c])+",")
            outfile1.write("\n")
    
    countryoutfile = open("Helpers/countries.txt", "w")
    for x in countries:
        countryoutfile.write(x+"\n")
    








