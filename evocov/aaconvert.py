import sys
import os
import numpy as np
import csv
import math
from tqdm import tqdm

def aaconvert(infile, outfile):
    reference = open("Data/spike_NT.txt", "r").readlines()
    reference = list(reference[0].upper())
    codondict = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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


    def NTmuttoAA(desc):
        orig = desc[0]
        new = desc[-2]
        pos = int(desc[1:-2])
        position = int(desc[1:-2])-1
        AApos = math.ceil(pos/3)
        if round(pos%3,2) == 1:
            codon = "".join(reference[position:(position+3)])
            newcodon = reference[position:(position+3)]
            newcodon[0] = new
            newcodon = "".join(newcodon)
        elif round(pos%3,2) == 2:
            codon = "".join(reference[position-1:position+2])
            newcodon = reference[position-1:position+2]
            newcodon[1] = new
            newcodon = "".join(newcodon)
        elif round(pos%3,2) == 0:
            codon = "".join(reference[position-2:position+1])
            newcodon = reference[position-2:position+1]
            newcodon[2] = new
            newcodon = "".join(newcodon)
        else:
            print(pos%3)

        try:
            oldAA = codondict[codon]
        except:
            oldAA = "X"
        try:
            newAA = codondict[newcodon]
        except:
            newAA = "X"

        if oldAA != newAA:
            final = oldAA+str(AApos)+newAA
            return final
        else:
            return "NA"





    NTdifffile = open("Data/"+infile, "r").readlines()
    AAdiffile = open("Data/"+outfile, "w")
    count = 0
    for line in tqdm(NTdifffile):
        if list(line)[0] == ">":
            count += 1
            AAdiffile.write(line)
            posns = []
        else:
            pos = math.ceil(int(line[1:-2])/3)
            if pos not in posns:
                AAmut = NTmuttoAA(line)
                if AAmut != "NA":
                    AAdiffile.write(AAmut+"\n")
                posns.append(pos)


    print(count)

