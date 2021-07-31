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
    for nt in range((AA*3)-3, (AA*3)):
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
    return e
    

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

def scoring(candidate, NTtable, AccNTtable, VarTables, TimTables, mutrate):
    codondict = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    AAscores = {
    "A":  5.0, "C":  7.0, "D":  10.0,
    "E": 10.0, "F": 5.0, "G":  5.0,
    "H": 10.0, "I": 5.0, "K": 10.0,
    "L": 5.0, "M": 5.0, "N":  8.0,
    "P":  7.0, "Q": 8.0, "R": 10.0,
    "S":  8.0, "T":  8.0, "V": 5.0,
    "W": 5.0, "Y": 5.0,
    }
    reference = open("Data/spikeAAcopy.txt", "r").readlines()
    reference = list(reference[0])
    NTreference = open("Data/spike_NT.txt", "r").readlines()
    NTreference = list(NTreference[0])
    indexs = (candidate.split(",_")[0].split(","))
    indexes = []
    AAs = []
    for i in indexs:
        indexes.append(int(i))
        AAs.append(reference[int(i)])
    if len(AAs)<5:
        return "NA"
    
    #Distance
    distance = float(candidate.split("_")[1].strip("\n"))
    distancescore = 5*((20-distance)/20)
    
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
    entro = []
    for x in indexes:
        entro.append(AAentropy(x, Last3mostable))
    entropies = mean(entro)
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
    entrop = mean(entropies)
    timescore = 5 - (5*entrop)

    #Flag multiple deletions
    deletions = []
    for x in indexes:
        if NTtable[x, -1] > 0.000005:
            deletions.append(x)
    
    
    totalscore = distancescore+lenscore+AAscore+locationscore+entropyscore+timescore+variantscore
    indexes = [x+1 for x in indexes]
    entropies = [str(x) for x in entro]
    if len(deletions) > 1:
        deletions = [str(x) for x in deletions]
        deletions = ",".join(deletions)


    #######Comparison of NT>NT mutation rates at each position to the baseline#########
    MutRateDict = {'A>C':0.039*mutrate, 'A>G':0.310*mutrate, 'A>T':0.123*mutrate,
                   'C>A':0.140*mutrate, 'C>G':0.022*mutrate, 'C>T':3.028*mutrate,
                   'G>A':0.747*mutrate, 'G>C':0.113*mutrate, 'G>T':2.953*mutrate,
                   'T>A':0.056*mutrate, 'T>C':0.261*mutrate, 'T>G':0.036*mutrate}
    MutRateDict2 = {'A>C':0.039, 'A>G':0.310, 'A>T':0.123,
                   'C>A':0.140, 'C>G':0.022, 'C>T':3.028,
                   'G>A':0.747, 'G>C':0.113, 'G>T':2.953,
                   'T>A':0.056, 'T>C':0.261, 'T>G':0.036}

    NTentropies = []
    NTposns = []
    NTs = []
    for x in indexes:
        result = entropy(x, AccNTtable)
        for y in result:
            NTentropies.append(y)
        for i in range((x*3)-3, (x*3)):
            NTposns.append(i)
            NTs.append(NTreference[i])
    bases = ["A","C","T","G"]
    predmut = False
    spentpos = []
    print("".join(AAs))
    print("".join(NTs))
    while predmut == False:
        maximumentropy = max(NTentropies)
        minentropyindex = NTentropies.index(maximumentropy)
        riskpos = NTposns[minentropyindex]
        riskAA = NTs[minentropyindex]
        spentNTs = []
        while len(spentNTs) != 3:
            highest = -math.inf
            for x in range(0,4):
                if bases[x] == riskAA:
                    continue
                if bases[x] in spentNTs:
                    continue
                else:
                    observed = AccNTtable[riskpos,x]
                    mutation = riskAA+">"+bases[x]
                    expected = MutRateDict[mutation]
                    if (observed - expected) > highest:
                        highest = observed-expected
                        sub = bases[x]
            print(riskpos)
            print(riskAA)
            print(round((riskpos % 3), 2))
            if round((riskpos % 3), 2) == 0:
                print("pos1")
                oldcodon = NTreference[riskpos]+NTreference[riskpos+1]+NTreference[riskpos+2]
                newcodon = sub+NTreference[riskpos+1]+NTreference[riskpos+2]
            elif round((riskpos % 3), 2) == 1:
                print("pos2")
                oldcodon = NTreference[riskpos-1]+NTreference[riskpos]+NTreference[riskpos+1]
                newcodon = NTreference[riskpos-1]+sub+NTreference[riskpos+1]
            elif round((riskpos % 3), 2) == 2:
                print("pos3")
                oldcodon = NTreference[riskpos-2]+NTreference[riskpos-1]+NTreference[riskpos]
                newcodon = NTreference[riskpos-2]+NTreference[riskpos-1]+sub
            
            if codondict[oldcodon] != codondict[newcodon]:
                predmut = True
                likelymutAA = codondict[oldcodon]+">"+str(round((riskpos-1)/3)+1)+">"+codondict[newcodon]
                likelymutNT = NTreference[riskpos]+">"+str(riskpos+1)+">"+sub
                spentNTs = ["A","C","T"]
            else:
                predmut = False
                spentNTs.append(sub)
        NTentropies.remove(maximumentropy)
        NTposns.remove(riskpos)
        NTs.remove(riskAA)
    
    if totalscore > 50:
        AAs = "".join(AAs)
        sitewiseentropies = ",".join(entropies)
        finallist = [AAs, indexes, distancescore, lenscore, AAscore, locationscore, entropyscore, timescore, variantscore, totalscore, sitewiseentropies, deletions, likelymutNT, likelymutAA]
        return finallist
    else:
        return "NA"

