import sys
import os
import numpy as np
import csv
import math
from evocov.makediffNT import seqparser
from evocov.aaconvert import aaconvert
from evocov.countersimple import simplecounter
from tqdm import tqdm

inputstring = input("Welcome to evocov, the SARS-CoV-2 epitope selection pipeline! First please pass the file paths of your latest GISAID FASTA file and the corresponding metadata, separated by a space.\n")
fasta = inputstring.split(" ")[0]
meta = inputstring.split(" ")[1]
correct = True
if not os.path.isfile(fasta) and os.path.isfile(meta):
    correct = False
while correct == False:
    inputstring = input("Uh oh! One of the file paths you entered was invalid. Please try again.\n")
    fasta = inputstring.split(" ")[0]
    meta = inputstring.split(" ")[1]
    if os.path.isfile(fasta) and os.path.isfile(meta):
        correct = True

metadata = list(csv.reader(open(meta, "r"), delimiter = "\t"))

NTfile = input("\n Next we will create a diff file to document the differences between each sequence and the SARS-CoV-2 reference, and parse the metadata of interest. Please type the desired name for this file.\n")
if os.path.isfile("Data/"+NTfile):
    repeat = input("\nThere is already a file with that name! Would you like to skip this step of the pipeline, completely rebuild this file from scratch, or update the file with any new sequences?(skip/scratch/update)\n")
    if repeat == "update":
        print("Cataloging the sequences already in "+NTfile)
        complete = []
        file = open("Data/"+NTfile, "r").readlines()
        for l in tqdm(file):
            if list(l)[0] == ">":
                acc = l.split("|")[0].split("> ")[1]
                complete.append(acc)
        complete = set(complete)
        seqparser(fasta, metadata, NTfile, complete)
    elif repeat == "scratch":
        complete = []
        seqparser(fasta, metadata, NTfile, complete)     
else:
    complete = []
    seqparser(fasta, metadata, NTfile, complete)


AAfile = input("\n Please enter the name desired for the Amino Acid version of the diff file.\n")
if os.path.isfile("Data/"+AAfile):
    repeat = input("\nThere is already a file with that name! Would you like to skip this step of the pipeline, or rebuild this file from scratch?(skip/scratch)\n")
    if repeat == "scratch":
        aaconvert(NTfile, AAfile)
else:
    aaconvert(NTfile, AAfile)


print("The pipeline will now proceed with counting and analysis and return representations of the mutational landscape and recommended epitope candidates")
simplecounter("NT", NTfile, "simplecountsNT.csv")
