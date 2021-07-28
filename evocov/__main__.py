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
from evocov.scoring import AAentropy
from evocov.scoringslidingwindow import scoringslidingwindow
from tqdm import tqdm
import json
import wget
import requests
import datetime
import subprocess
from itertools import repeat
import operator
from operator import itemgetter
from evocov.epitopefrequencycalculator import freqcalc

backdoor = False
try:
    if sys.argv[1] == "backdoor":
        backdoor = True
        default = True
    else:
        backdoor = False
        default = True
except:
    backdoor = False
    default = False
    
if backdoor == False:
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
    if not (os.path.isfile(fasta) and os.path.isfile(meta)):
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
        if backdoor == True:
            repeat = "skip"
            print("\nWe have detected an existing NT diff file, so this step will be skipped.\n")
        else:
            repeat = "scratch"
            print("\nWe have detected an existing NT diff file, so existing sequences in the diff file will be catalogued and only new ones will be added.\n")
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
        if backdoor == True:
            repeat == "skip"
        else:
            repeat = "scratch"
            print("\nMaking an Amino Acid Version of the Diff File!\n")
    else:
        repeat = input("\nThere is already a file with that name! Would you like to skip this step of the pipeline, or rebuild this file from scratch?(skip/scratch)\n")
    if repeat == "scratch":
        aaconvert(NTfile, AAfile)
else:
    print("\nMaking an Amino Acid Version of the Diff File!\n")
    aaconvert(NTfile, AAfile)

#Counting
subprocess.call("rm -r Analysis/date", shell = True)
subprocess.call("rm -r Analysis/variant", shell = True)
if not os.path.isdir("Analysis"):
    os.makedirs("Analysis")
print("The pipeline will now proceed with counting and analysis of the mutational landscape of the SARS-CoV-2 Genome.")

if backdoor == True:
    numsequences = 100
    NTcounting = simplecounter("NT", NTfile, "Analysis/simplecountsNT.csv")
    numsequences = NTcounting[0]
    mutrate = NTcounting[1]
else:
    print("Round 1/4")
    NTcounting = simplecounter("NT", NTfile, "Analysis/simplecountsNT.csv")
    numsequences = NTcounting[0]
    mutrate = NTcounting[1]
if backdoor == True:
    NTtable = np.genfromtxt('Analysis/simplecountsAA.csv', delimiter=',')
else:
    print("Round 2/4")
    NTtable = simplecounter("AA", AAfile, "Analysis/simplecountsAA.csv")
print("Round 3/4")
timtables = metasplitcounter("AA", AAfile, "date", ["2019-12", "2020-01","2020-02","2020-03","2020-04","2020-05","2020-06","2020-07","2020-08","2020-09","2020-10","2020-11","2020-12","2021-01","2021-02","2021-03","2021-04","2021-05","2021-06"])


if default == True:
    print("Because you are running on default, the variants used for this analysis will be those classified as VOC: B.1.1.7, B.1.351, B.1.427, B.1.429, P.1, B.1.617.2")
    varguments = ["B.1.1.7", "B.1.351", "B.1.427", "B.1.429","P.1", "B.1.617.2",  "C.37"]
else:
    var = input("Please enter the names of the variants you wish to use in this analysis separated by spaces. If you wish to proceed with the default variants (B.1.1.7, B.1.351, B.1.427, B.1.429, P.1, B.1.617.2), just press enter")
    if input == "":
        varguments = ["B.1.1.7", "B.1.351", "B.1.427", "B.1.429","P.1", "B.1.617.2", "C.37"]
    else:
        varguments = var.split(" ")
print("Round 4/4")        
vartables = metasplitcounter("AA", AAfile, "variant", varguments)

AccNTtable = np.genfromtxt('Analysis/simplecountsNT.csv', delimiter=',')

#Scoring
print("\nLoading in and scoring epitope candidates\n")

candidates = open("Data/vendruscolocandidates.txt").readlines()
scores = []
for x in tqdm(candidates):
    score = scoring(x, NTtable, AccNTtable, vartables, timtables, mutrate)
    if score != "NA":
        scores.append(score)

sortedscores = sorted(scores, key = itemgetter(9), reverse = True)

freqcalc(NTfile, candidates, ['2021-07'])

with open("Analysis/scoredepitopes.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(sortedscores)


'''         
candidates = open("Data/slidingwindowcandidates.txt").readlines()
scores = []
for x in tqdm(candidates):
    score = scoringslidingwindow(x, NTtable, vartables, timtables)
    if score != "NA":
        scores.append(score)

sortedscores = sorted(scores, key = itemgetter(9), reverse = True)
with open("Analysis/scoredslidingwindowepitopes.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(sortedscores)
'''   

#Download current case data for all countries and normalise
print("Retrieving current case data for all countries from the WHO, this will be used to normalise the counts by country.\n")
wget.download("https://covid19.who.int/WHO-COVID-19-global-table-data.csv", "WHOcasedata.csv")


#Make PDF
if default == True:
    if not os.path.isdir("Plots"):
        os.makedirs("Plots")
    print("\nCreating a PDF of the pipeline results")
    subprocess.call("Rscript evocov/plotter.R Analysis/simplecountsNT.csv "+str(numsequences)+" Analysis/simplecountsAA.csv", shell = True)

else:
    plotornot = input("Would you like to use R to pipe the output to a PDF? (y/n)")
    if plotornot == "y":
        if not os.path.isdir("Plots"):
            os.makedirs("Plots")
        print("\nCreating a PDF of the pipeline results")
        subprocess.call("Rscript evocov/plotter.R Analysis/simplecountsNT.csv "+str(numsequences)+" Analysis/simplecountsAA.csv", shell = True)
        

subprocess.call("rm WHO*", shell = True)

try:
    subprocess.call("curl http://textbelt.com/text -d number="+sys.argv[3]+" -d message=\"Pipeline Analysis Complete!\" -d key=e237256cdcb7af72df888a7558f92c0e97b0fb55OVbGpYEsvYujZXDCVi0Rtvom6 ", shell = "True")
    sys.exit()
except:
    sys.exit()
