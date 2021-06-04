import os
import sys
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add
import numpy as np


filepath = "Data/2021-05-26_unmasked.fa"
reference = open("Data/spike_NT.txt", "r").readlines()
reference = list(reference[0].lower())
accessions = open("accessions.txt", "r").read().splitlines()
dates = open("quarters.txt", "r").read().splitlines()
strains = open("strains.txt", "r").read().splitlines()
locations = open("countries.txt", "r").read().splitlines()

accessions = [i.replace('"', '') for i in accessions]
strains = [i.replace('"', '') for i in strains]
locations = [i.replace('"', '') for i in locations]

def mutationchecker(num):
    if reference[num] == sequence[num]:
        return False
    else:
        if sequence[num] != "n":
            return True
        else:
            return False


bases = ["a","c","t","g","n","-"]

with open(filepath, mode = "r") as handle:
    outfile = open("difffile.txt", "w")
    errorfile = open("misaligned.txt", "w")
    count = 0
    errors = 0
    for record in tqdm(SeqIO.parse(handle, 'fasta')):
        accession = record.id
        try:
            index = accessions.index(str(accession))
        except ValueError:
            continue
        time = int(dates[index])
        variant = strains[index]
        location = locations[index]
        sequence = list(record.seq)
        sequence = sequence[21562:25384]
        try:
            mutations = list(filter(mutationchecker, range(0, (len(reference)-1))))
        except:
            continue
        
        label = "> "+accession+"|"+str(time)+"|"+variant+"|"+location

        if len(mutations) > 30:
            errors += 1
            errorfile.write(label+"\n")
            for x in sequence:
                errorfile.write(x)
            errorfile.write("\n")
        elif len(mutations) > 0:
            count += 1
            outfile.write(label+"\n")
            for x in mutations:
                outfile.write(reference[x].upper()+str(x)+sequence[x].upper()+"\n")
        else:
            count += 1
            outfile.write(label+"\n")
            
        '''        
        if count > 100:
            outfile.close()
            errorfile.close()
            sys.exit()
        '''    

print(errors)
print(count)
       
outfile.close()
errorfile.close()

subprocess.call("curl http://textbelt.com/text -d number=353877910680 -d message=\"Analysis Complete, "+str(count)+" Sequences Analysed\" -d key=e237256cdcb7af72df888a7558f92c0e97b0fb55OVbGpYEsvYujZXDCVi0Rtvom6 ", shell = "True")

