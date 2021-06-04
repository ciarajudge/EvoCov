import os
import sys
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add
import numpy as np


filepath = "2021-05-26_unmasked.fa"
reference = open("spike_nt.txt", "r").readlines()
reference = list(reference[0].lower())
score = [0]*len(reference)


bases = ["a","c","t","g","n","-"]

accessions = open("accessions.txt", "r").read().splitlines()
dates = open("quarters.txt", "r").read().splitlines()
strains = open("strains.txt", "r").read().splitlines()
accessions = [i.replace('"', '') for i in accessions] 
print(accessions[1:20])
print(accessions[1])

#Not the most logical way to do this but whatever
_2019Q4matrix = np.zeros((len(reference), 6))
_2020Q1matrix = np.zeros((len(reference), 6))
_2020Q2matrix = np.zeros((len(reference), 6))
_2020Q3matrix = np.zeros((len(reference), 6))
_2020Q4matrix = np.zeros((len(reference), 6))
_2021Q1matrix = np.zeros((len(reference), 6))
_2021Q2matrix = np.zeros((len(reference), 6))
_2019Q4 = 0
_2020Q1 = 0
_2020Q2 = 0
_2020Q3 = 0
_2020Q4 = 0
_2021Q1 = 0
_2021Q2 = 0
#Variant tables
B117matrix = np.zeros((len(reference), 6))
B1135matrix = np.zeros((len(reference), 6))
B1427matrix = np.zeros((len(reference), 6))
B1429matrix = np.zeros((len(reference), 6))
P1matrix = np.zeros((len(reference), 6))
B1617matrix = np.zeros((len(reference), 6))
B117 = 0
B1135 = 0
B1427 = 0
B1429 = 0
P1 = 0
B1617 = 0

variants = ["B.1.1.7", "B.1.315", "B.1.427", "B.1.429", "P.1", "B.1.617.2"]
times = [20194, 20201, 20202, 20203, 20204, 20211, 20212]
try:
    with open(filepath, mode = "r") as handle:
        count = 0
        for record in tqdm(SeqIO.parse(handle, 'fasta')):
            accession = record.id
            try:
                index = accessions.index(str(accession))
            except ValueError:
                continue
    
            if dates[index] == 'NA':
                continue
            if strains[index] == 'NA':
                continue
            time = int(dates[index])
            variant = strains[index]
            if time == 20200:
                continue
            test = time in times or variant in variants
            if test == False:
                continue
            count += 1
            tempmatrix = np.zeros((len(reference), 6))

            sequence = list(record.seq)
            sequence = sequence[21562:25384]
        
            for nt in range((len(sequence)-1)):
                if sequence[nt] != reference[nt]:
                    if sequence[nt] in bases:
                        row = nt
                        col = bases.index(sequence[nt])
                        tempmatrix[row,col] = tempmatrix[row,col] + 1
                    else:
                        continue

            
            #Update for time
            if time == 20194:
                _2019Q4matrix = np.add(_2019Q4matrix, tempmatrix)
                _2019Q4 += 1
            elif time == 20201:
                _2020Q1matrix = np.add(_2020Q1matrix, tempmatrix)
                _2020Q1 += 1
            elif time == 20202:
                _2020Q2matrix = np.add(_2020Q2matrix, tempmatrix)
                _2020Q2 += 1
            elif time == 20203:
                _2020Q3matrix = np.add(_2020Q3matrix, tempmatrix)
                _2020Q3 += 1
            elif time == 20204:
                _2020Q4matrix = np.add(_2020Q4matrix, tempmatrix)
                _2020Q4 += 1
            elif time == 20211:
                _2021Q1matrix = np.add(_2021Q1matrix, tempmatrix)
                _2021Q1 += 1
            else:
                _2021Q2matrix = np.add(_2021Q2matrix, tempmatrix)
                _2021Q2 += 1
            
           #Update for variant
            if variant == "B.1.1.7":
                B117matrix = np.add(B117matrix, tempmatrix)
                B117 += 1
            elif variant == "B.1.315":
                B1135matrix = np.add(B1135matrix, tempmatrix)
                B1135 += 1
            elif variant == "B.1.427":
                B1427matrix = np.add(B1427matrix, tempmatrix)
                B1427 += 1
            elif variant == "B.1.429":
                B1429matrix = np.add(B1429matrix, tempmatrix)
                B1429 += 1
            elif variant == "P.1":
                P1matrix = np.add(P1matrix, tempmatrix)
                P1 += 1
            elif variant == "B.1.617.2":
                B1617matrix = np.add(B1617matrix, tempmatrix)
                B1617 += 1

            if count > 20000:
                break

except:
    subprocess.call("curl http://textbelt.com/text -d number=353877910680 -d message=\"Something's gone wrong :(\" -d key=e237256cdcb7af72df888a7558f92c0e97b0fb55OVbGpYEsvYujZXDCVi0Rtvom6 ", shell = "True")




#Convert matrices to frequencies
'''
_2019Q4matrix = _2019Q4matrix/_2019Q4
_2020Q1matrix = _2020Q1matrix/_2020Q1
_2020Q2matrix = _2020Q2matrix/_2020Q2
_2020Q3matrix = _2020Q3matrix/_2020Q3
_2020Q4matrix = _2020Q4matrix/_2020Q4
_2021Q1matrix = _2021Q1matrix/_2021Q1
_2021Q2matrix = _2021Q1matrix/_2021Q2
'''



np.savetxt("2019Q4.csv", _2019Q4matrix, delimiter=",")
np.savetxt("2020Q1.csv", _2020Q1matrix, delimiter=",")
np.savetxt("2020Q2.csv", _2020Q2matrix, delimiter=",")
np.savetxt("2020Q3.csv", _2020Q3matrix, delimiter=",")
np.savetxt("2020Q4.csv", _2020Q4matrix, delimiter=",")
np.savetxt("2021Q1.csv", _2021Q1matrix, delimiter=",")
np.savetxt("2021Q2.csv", _2021Q2matrix, delimiter=",")

np.savetxt("B117.csv", B117matrix, delimiter=",")
np.savetxt("B1135.csv", B1135matrix, delimiter=",")
np.savetxt("B1427.csv", B1427matrix, delimiter=",")
np.savetxt("B1429.csv", B1429matrix, delimiter=",")
np.savetxt("P1.csv", P1matrix, delimiter=",")
np.savetxt("B1617.csv", B1617matrix, delimiter=",")

subprocess.call("curl http://textbelt.com/text -d number=353877910680 -d message=\"Analysis Complete, "+str(count)+" Sequences Analysed\" -d key=e237256cdcb7af72df888a7558f92c0e97b0fb55OVbGpYEsvYujZXDCVi0Rtvom6 ", shell = "True")
print("years")
print(str(_2019Q4))
print(str(_2020Q1))
print(str(_2020Q2))
print(str(_2020Q3))
print(str(_2020Q4))
print(str(_2021Q1))
print(str(_2021Q2))

print("strains")
print(str(B117))
print(str(B1135))
print(str(B1427))
print(str(B1429))
print(str(P1))
print(str(B1617))




