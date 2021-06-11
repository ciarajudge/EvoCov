import os
import sys
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add
import numpy as np
import csv
from tabulate import tabulate


filepath = "Data/2021-06-11_unmasked.fa"
reference = open("Data/spike_NT.txt", "r").readlines()
reference = list(reference[0].lower())
metadata = list(csv.reader(open("Data/metadata.tsv", "r"), delimiter = "\t"))

dictionary = {item[2]:item[3:5]+item[10:12] for item in metadata}

def mutationchecker(num):
    if reference[num] == sequence[num]:
        return False
    else:
        return True

bases = ["a","c","t","g","n","-"]

with open(filepath, mode = "r") as handle:
    outfile = open("difffile2.txt", "w")
    errorfile = open("misaligned2.txt", "w")
    count = 0
    errors = 0
    for record in tqdm(SeqIO.parse(handle, 'fasta')):
        accession = record.id
        try:
            meta = dictionary[accession]
        except ValueError:
            continue
        sequence = list(record.seq)
        sequence = sequence[21562:25384]
        try:
            mutations = list(filter(mutationchecker, range(0, (len(reference)-1))))
        except:
            continue
        
        label = "> "+accession+"|"+str(meta[0])+"|"+meta[3]+"|"+meta[2]+"|"+meta[1]

        if len(mutations) > 500:
            errors += 1
            errorfile.write(label+"\n")
            for x in sequence:
                errorfile.write(x)
            errorfile.write("\n")
        elif len(mutations) > 0:
            count += 1
            outfile.write(label+"\n")
            for x in mutations:
                outfile.write(reference[x].upper()+str(x+1)+sequence[x].upper()+"\n")
        else:
            count += 1
            outfile.write(label+"\n")
            
        """        
        if count > 100:
            outfile.close()
            errorfile.close()
            sys.exit()
        """   

print(errors)
print(count)
       
outfile.close()
errorfile.close()

subprocess.call("curl http://textbelt.com/text -d number=353877910680 -d message=\"Analysis Complete, "+str(count)+" Sequences Analysed\" -d key=e237256cdcb7af72df888a7558f92c0e97b0fb55OVbGpYEsvYujZXDCVi0Rtvom6 ", shell = "True")

