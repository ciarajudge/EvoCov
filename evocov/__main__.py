import sys
import os
import numpy as np
import csv
import math
from evocov.makediffNT import seqparser
from evocov.aaconvert import aaconvert
from tqdm import tqdm

inputstring = input("Welcome to evocov, the SARS-CoV-2 epitope selection pipeline! First please pass the filepaths of your latest GISAID FASTA file and the corresponding metadata, separated by a space.\n")
fasta = inputstring.split(" ")[0]
meta = inputstring.split(" ")[1]
metadata = list(csv.reader(open(meta, "r"), delimiter = "\t"))

outfile = input("\n Next we will create a diff file to document the differences between each sequence and the SARS-CoV-2 reference, and parse the metadata of interest. Please type the desired name for this file. If another file of the same name exists, it will be checked for already analysed accession numbers, and these will not be repeated. Analysis from scratch takes approximately 2 hours, and a weekly update will take approximately 10 minutes.\n")

complete = []
if os.path.isfile("Data/"+outfile):
    file = open("Data/"+outfile, "r").readlines()
    for l in tqdm(file):
        if list(l)[0] == ">":
            acc = l.split("|")[0].split("> ")[1]
            complete.append(acc)
complete = set(complete)


print("\n Already completed sequences catalogued, GISAID Fasta will now be parsed.")
seqparser(fasta, metadata, outfile, complete)

AAfile = input("\n Please enter the name desired for the Amino Acid version of the diff file.\n")
aaconvert(outfile, AAfile)

