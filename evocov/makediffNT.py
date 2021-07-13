import os
import sys
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add


def seqparser(filepath, metadata, otfile, complete):
    def mutationchecker(num):
        if reference[num] == sequence[num]:
            return False
        else:
            return True
    reference = open("Data/spike_NT.txt", "r").readlines()
    reference = list(reference[0].lower())
    dictionary = {item[2]:item[3:5]+item[10:12] for item in metadata}
    with open(filepath, mode = "r") as handle:
        if len(complete) < 1:
            outfile = open("Data/"+otfile, "w")
        else:
            outfile = open("Data/"+otfile, "a")
        errorfile = open("Data/misaligned.txt", "w")
        bases = ["a","c","t","g","n","-"]
        count = 0
        errors = 0
        added = 0
        for record in tqdm(SeqIO.parse(handle, 'fasta')):
            accession = record.id
            if accession in complete:
                count += 1
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
                print("error")
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

            added += 1  

    print(str(errors)+" sequences with high diff count detected and saved to misaligned.txt")
    print(str(count)+" sequences added to "+otfile)
       
    outfile.close()
    errorfile.close()
    return count

    #subprocess.call("curl http://textbelt.com/text -d number=353877910680 -d message=\"Analysis Complete, "+str(added)+" Sequences Added\" -d key=e237256cdcb7af72df888a7558f92c0e97b0fb55OVbGpYEsvYujZXDCVi0Rtvom6 ", shell = "True")

