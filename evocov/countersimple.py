import sys
import os
import numpy as np
import csv
import math
from tqdm import tqdm

def simplecounter(form, infile, outfile):
    infile = open("Data/"+infile, "r").readlines()
    if form == "NT" :
        reference = open("Data/spike_NT.txt", "r").readlines()
        reference = list(reference[0].lower())
        bases = ["A","C","T","G","N","-"]
        posmutmatrix = np.zeros((len(reference), len(bases)))
    elif form == "AA" :
        reference = open("Data/spike_AA.txt", "r").readlines()
        reference = list(reference[0].lower())
        bases = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P","S", "T", "W", "Y", "V", "X", "*", "_"]
        posmutmatrix = np.zeros((len(reference)+10, len(bases)))

    count = 0
    otherbases = []
    for line in tqdm(infile):
        if list(line)[0] == ">":
            count += 1
            posns = []
        else:
            base = line[-2]
            position = int(line[1:-2])
            if position not in posns:
                if base in bases:
                    row = position
                    col = bases.index(base)
                    posmutmatrix[row,col] = posmutmatrix[row,col] + 1
                    posns.append(position)
                else:
                    if base in otherbases:
                        continue
                    else:
                        otherbases.append(base)

    print(count)
    print(otherbases)
    np.savetxt(outfile, (posmutmatrix/count), delimiter=",")
    return(posmutmatrix/count)

