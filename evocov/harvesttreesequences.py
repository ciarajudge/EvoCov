import os
import sys
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add


accessions = open("subsetlabels.txt", "r").readlines()
accessions = [x.strip("\n") for x in accessions]
print(accessions)
outfile = open("sequences.nuc", "a")
filepath = "2021-06-25_unmasked.fa"
with open(filepath, mode = "r") as handle:
    for record in tqdm(SeqIO.parse(handle, 'fasta')):
        accession = record.id
        if accession in set(accessions):
            outfile.write(accession+"\n")
            sequence = list(record.seq)[21562:25384]
            sequence = "".join(sequence)
            outfile.write(sequence+"\n")

outfile.close()



