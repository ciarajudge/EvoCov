import sys
import os
import numpy as np
import csv
import math
from tqdm import tqdm

def mutationmaker(form, loci, newbases, outfilename1, outfilename2):
    if form == "NT" :
        reference = open("Data/spike_NT.txt", "r").readlines()
        reference = list(reference[0].lower())
    elif form == "AA" :
        reference = open("Data/spike_AA.txt", "r").readlines()
        reference = list(reference[0].lower())

    newseq = reference
    for x in range(0, len(loci)):
        newseq[loci[x]-1] = newbases[x]

    newseqRBD = newseq[318:541]
    outseq = "".join(newseq).upper()
    outseqRBD = "".join(newseqRBD).upper()
    
    outfile1 = open(outfilename1, "w")
    outfile1.write(outseq)

    outfile2 = open(outfilename2, "w")
    outfile2.write(outseqRBD)
        

mutationmaker("AA", [493], ["E"], "spikeAA_Q493E.fa", "RBDAA_Q493E.fa")
mutationmaker("AA", [452], ["V"], "spikeAA_L452V.fa", "RBDAA_L452V.fa")
