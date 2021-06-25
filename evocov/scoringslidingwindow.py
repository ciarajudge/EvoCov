import sys
import os
import numpy as np
import csv
import math
from tqdm import tqdm
from statistics import mean
from statistics import variance


def entropy(AA, NTtable):
    e = []
    for nt in range(AA*3, (AA*3)+4):
        total = NTtable[nt,0]+NTtable[nt,1]+NTtable[nt,2]+NTtable[nt,3]+NTtable[nt,5]
        pA = NTtable[nt,0]/total
        pC = NTtable[nt,1]/total
        pT = NTtable[nt,2]/total
        pG = NTtable[nt,3]/total
        pGap = NTtable[nt,5]/total
        try:
            A = pA*math.log(pA)
        except ValueError:
            A = 0
        try:
            C = pC*math.log(pC)
        except ValueError:
            C = 0
        try:
            T = pT*math.log(pT)
        except ValueError:
            T = 0
        try:
            G = pG*math.log(pG)
        except ValueError:
            G = 0
        try:
            Gap = pGap*math.log(pGap)
        except ValueError:
            Gap = 0
        entropy = -(A+C+T+G+Gap)
        e.append(entropy/math.log(5))
    codone = mean(e)
    return codone
    

def AAentropy(AA, NTtable):
    ps = []
    total = (sum(NTtable[AA,]) - NTtable[AA, 20]) - NTtable[AA, 21]
    ranges = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22]
    for base in ranges:
        p = NTtable[AA,base]/total
        try:
            plog = p*math.log(p)
        except:
            plog = 0
        ps.append(plog)
    entropy = -(sum(ps))
    normentropy = entropy/math.log(21)
    return normentropy

def scoringslidingwindow(candidate, NTtable, VarTables, TimTables):
    externals = open("Data/externals.txt", "r").readlines()
    AAscores = {
    "A":  5.0, "C":  7.0, "D":  10.0,
    "E": 10.0, "F": 5.0, "G":  5.0,
    "H": 10.0, "I": 5.0, "K": 10.0,
    "L": 5.0, "M": 5.0, "N":  8.0,
    "P":  7.0, "Q": 8.0, "R": 10.0,
    "S":  8.0, "T":  8.0, "V": 5.0,
    "W": 5.0, "Y": 5.0,
    }
    reference = open("Data/spike_AA.txt", "r").readlines()
    reference = list(reference[0])
    indexs = (candidate.split(",_")[0].split(","))
    indexes = []
    AAs = []
    for i in indexs:
        indexes.append(int(i))
        AAs.append(reference[int(i)])
    if len(AAs)<5:
        return "NA"
    
    #Internal
    external = 0
    for i in indexes:
        if i in set(externals):
            external += 1
    external = external/7
    externalscore = external*5
    
    #Length
    if len(indexes)>5:
        lenscore = 5
    else:
        lenscore = len(indexes)
       
    #AA Score
    AAscore = []
    for x in AAs:
        AAscore.append(AAscores[x])
    AAscore = (mean(AAscore)/10)*10
    
    #Location
    inRBD = 0
    for x in indexes:
        if x in range(318, 541):
            inRBD += 1
    inRBD = inRBD/len(indexes)
    RBDscore = 5*inRBD
    locationscore = RBDscore
    
    #Entropy
    Last3mostable = TimTables[:,:,-1]+TimTables[:,:,-2]+TimTables[:,:,-3]
    entropies = []
    for x in indexes:
        entropies.append(AAentropy(x, Last3mostable))
    entropies = mean(entropies)
    entropyscore = 65 - (65*entropies)
    
    
    #Uniform Entropy across Variants
    entropyvar = []
    for x in indexes:
        entropies = []
        for v in range(0, (np.size(VarTables, 2))):    
            ent = AAentropy(x, VarTables[:,:,v])
            entropies.append(ent)
        entropyvar.append(variance(entropies))
    entropyvariance = mean(entropyvar)
    variantscore = (1 - (entropyvariance/0.7))*5
  
    #Uniform Entropy across Time
    entropies = []
    for x in indexes:
        entropies.append(AAentropy(x,NTtable))
    entropies = mean(entropies)
    timescore = 5 - (5*entropies)
    
    
    totalscore = externalscore+lenscore+AAscore+locationscore+entropyscore+timescore+variantscore
    indexes = [x+1 for x in indexes]
    if totalscore > 50:
        AAs = "".join(AAs)
        finallist = [AAs, indexes, externalscore, lenscore, AAscore, locationscore, entropyscore, timescore, variantscore, totalscore]
        return finallist
    else:
        return "NA"

