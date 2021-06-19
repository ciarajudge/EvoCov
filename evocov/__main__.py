import sys
import os
import numpy as np
import csv
import math
from evocov.makediffNT import seqparser
from evocov.aaconvert import aaconvert
from evocov.countersimple import simplecounter
from evocov.metasplitcounter import metasplitcounter
from evocov.scoring import scoring
from evocov.scoring import entropy
from tqdm import tqdm
import json
import wget
import requests
import pdfkit
import datetime
import subprocess
from itertools import repeat
import operator
from operator import itemgetter

#Initialise
if len(sys.argv) > 1:
    print("Welcome to evocov, the SARS-CoV-2 epitope selection pipeline! Because you have already passed your metadata and fasta file, the pipeline will run on default.")
    default = True
    fasta = sys.argv[1]
    meta = sys.argv[2]
else:
    default = False
    inputstring = input("Welcome to evocov, the SARS-CoV-2 epitope selection pipeline! First please pass the file paths of your latest GISAID FASTA file and the corresponding metadata, separated by a space.\n")
    fasta = inputstring.split(" ")[0]
    meta = inputstring.split(" ")[1]
    
#Make sure correct fasta and metadata files have been passed
correct = True
if not os.path.isfile(fasta) and os.path.isfile(meta):
    correct = False
while correct == False:
    inputstring = input("Uh oh! One of the file paths you entered was invalid. Please try again (enter the paths separated by a space).\n")
    fasta = inputstring.split(" ")[0]
    meta = inputstring.split(" ")[1]
    if os.path.isfile(fasta) and os.path.isfile(meta):
        correct = True
metadata = list(csv.reader(open(meta, "r"), delimiter = "\t"))

#Get NT diff file or file path to existing one
if default == True:
    NTfile = "difffile.txt"
else:
    NTfile = input("\nNext we will create a diff file to document the differences between each sequence and the SARS-CoV-2 reference, and parse the metadata of interest. Please type the desired name for this file.\n")
if os.path.isfile("Data/"+NTfile):
    if default == True:
        repeat = "skip"
        print("\nWe have detected an existing diff file, so this step will be skipped!\n")
    else:
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
        numsequences = seqparser(fasta, metadata, NTfile, complete)
    elif repeat == "scratch":
        complete = []
        numsequences = seqparser(fasta, metadata, NTfile, complete)        
else:
    complete = []
    numsequences = seqparser(fasta, metadata, NTfile, complete)

#Convert to AA diff file or get file path to existing one
if default == True:
    AAfile = "AAdifffile.txt"
else:
    AAfile = input("\n Please enter the name desired for the Amino Acid version of the diff file.\n")
if os.path.isfile("Data/"+AAfile):
    if default == True:
        repeat = "skip"
        print("\nWe have detected an existing AA diff file, so this step will be skipped!\n")
    else:
        repeat = input("\nThere is already a file with that name! Would you like to skip this step of the pipeline, or rebuild this file from scratch?(skip/scratch)\n")
    if repeat == "scratch":
        aaconvert(NTfile, AAfile)
else:
    aaconvert(NTfile, AAfile)

#Counting
if not os.path.isdir("Analysis"):
    os.makedirs("Analysis")
print("The pipeline will now proceed with counting and analysis and return representations of the mutational landscape and recommended epitope candidates.")
NTtable = simplecounter("NT", NTfile, "Analysis/simplecountsNT.csv")
#simplecounter("AA", AAfile, "Analysis/simplecountsAA.csv")
#varcounts = metasplitcounter("AA", AAfile, "variant", ["B.1.1.7", "B.1.351", "B.1.427", "B.1.429","P.1", "B.1.617.2"])


#Scoring
candidates = open("Data/candidates.txt").readlines()
scores = []
for x in tqdm(candidates):
    score = scoring(x, NTtable)
    if score != "NA":
        scores.append(score)
with open("scoredepitopes.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(scores)

sortedscores = sorted(scores, key = itemgetter(7))

with open("sortedscoredepitopes.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(sortedscores)
    
'''
#Download current case data for all countries and normalise
print("Retrieving current case data for all countries from the WHO, this will be used to normalise the counts by country.\n")
wget.download("https://covid19.who.int/WHO-COVID-19-global-table-data.csv", ".")

#Make PDF
if not os.path.isdir("Plots"):
    os.makedirs("Plots")
print("Creating a PDF of the pipeline results")
subprocess.call("Rscript evocov/plotter.R Analysis/simplecountsNT.csv "+str(numsequences)+" Analysis/simplecountsAA.csv variant", shell = True)

subprocess.call("rm WHO*", shell = True)
'''

























