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
        
mutationmaker("AA", [354], ["K"], "spikeAA_N354K.fa", "RBDAA_N354K.fa")
mutationmaker("AA", [354], ["S"], "spikeAA_N354S.fa", "RBDAA_N354S.fa")
mutationmaker("AA", [354], ["I"], "spikeAA_N354I.fa", "RBDAA_N354I.fa")
mutationmaker("AA", [500], ["P"], "spikeAA_T500P.fa", "RBDAA_T500P.fa")
mutationmaker("AA", [500], ["C"], "spikeAA_Y449C.fa", "RBDAA_Y449C.fa")
mutationmaker("AA", [500], ["E"], "spikeAA_Q493E.fa", "RBDAA_Q493E.fa")
