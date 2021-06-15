import os
import sys
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add
import numpy as np
import csv



filepath = sys.argv[1]
reference = open("Data/spike_NT.txt", "r").readlines()
reference = list(reference[0].lower())
metadata = list(csv.reader(open(sys.argv[2], "r"), delimiter = "\t"))
if os.path.isfile("Helpers/finished.txt"):
    complete = open("Helpers/finished.txt", "r").readlines()
    complete = set([x.split("\n")[0] for x in complete])
    completefile = open("Helpers/finished.txt", "a")
else:
    complete = []
    completefile = open("Helpers/finished.txt", "w")

outfile = open("Data/difffile", "a")
errorfile = open("Data/misaligned.txt", "w")


dictionary = {item[2]:item[3:5]+item[10:12] for item in metadata}

def mutationchecker(num):
    if reference[num] == sequence[num]:
        return False
    else:
        return True

bases = ["a","c","t","g","n","-"]

with open(filepath, mode = "r") as handle:
    count = 0
    errors = 0
    added = 0
    for record in tqdm(SeqIO.parse(handle, 'fasta')):
        accession = record.id
        if accession in complete:
            continue       
        try:
            meta = dictionary[accession]
        except:
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

        completefile.write(accession+"\n")
        added += 1
        """        
        if count > 100:
            outfile.close()
            errorfile.close()
            sys.exit()
        """   

print(errors)
print(count)
print(added)
       
outfile.close()
errorfile.close()
completefile.close()

subprocess.call("curl http://textbelt.com/text -d number=353877910680 -d message=\"Analysis Complete, "+str(added)+" Sequences Added\" -d key=e237256cdcb7af72df888a7558f92c0e97b0fb55OVbGpYEsvYujZXDCVi0Rtvom6 ", shell = "True")

