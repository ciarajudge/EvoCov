import sys
import os
import numpy as np
import csv
import math
from tqdm import tqdm
from statistics import mean


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
    

def scoring(candidate, NTtable):
    AAscores = {
    "A":  67.0, "C":  86.0, "D":  91.0,
    "E": 109.0, "F": 135.0, "G":  48.0,
    "H": 118.0, "I": 124.0, "K": 135.0,
    "L": 124.0, "M": 124.0, "N":  96.0,
    "P":  90.0, "Q": 114.0, "R": 148.0,
    "S":  73.0, "T":  93.0, "V": 105.0,
    "W": 163.0, "Y": 141.0,
    }
    reference = open("Data/spike_AA.txt", "r").readlines()
    reference = list(reference[0])
    indexs = (candidate.split(",_")[0].split(","))
    indexes = []
    AAs = []
    for i in indexs:
        indexes.append(int(i))
        AAs.append(reference[int(i)])
    
    #Distance
    distance = float(candidate.split("_")[1].strip("\n"))
    distancescore = 10*((15-distance)/15)
    #print(distancescore)
    #Length
    if len(indexes)>10:
        lenscore = 10
    else:
        lenscore = len(indexes)
    #print(lenscore)
    #AA Score
    AAscore = []
    for x in AAs:
        AAscore.append(AAscores[x])
    AAscore = (mean(AAscore)/170)*10
    #print(AAscore)
    #RBD
    inRBD = 0
    for x in indexes:
        if x in range(318, 541):
            inRBD += 1
    inRBD = inRBD/len(indexes)
    RBDscore = 10*inRBD
    #print(RBDscore)
    
    #Entropy
    entropies = []
    for x in indexes:
        entropies.append(entropy(x, NTtable))
    entropies = mean(entropies)
    entropyscore = 60 - (60*entropies)
    #print(entropyscore)          

    totalscore = distancescore+lenscore+AAscore+RBDscore+entropyscore
    #print(totalscore)
