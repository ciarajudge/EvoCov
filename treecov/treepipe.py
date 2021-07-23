import os
import sys
import re
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add


subprocess.call("Rscript subsampletree.R", shell = True)

files = os.listdir("labels")
files.remove(".DS_Store")

accessions = []
allaccessions = []
for f in files:
    accs = open("labels/"+f, "r").readlines()
    accessions.append([a.strip("\n") for a in accs])
    for x in accs:
        allaccessions.append(x)
allaccessions = [x.strip("\n") for x in allaccessions]

outfiles = []
outfilenames = []
for f in range(0, len(files)):
    outfiles.append(open("sequences/seqs"+str(f+1)+".nuc", "w"))
    outfilenames.append("sequences/seqs"+str(f+1)+".nuc")
    outfiles[f].write(" 1000 3822\n\n")

counts = [0,0,0,0,0,0,0,0,0,0]
totalcounts = 0
coveredacc = []
filepath = "2021-07-21_unmasked.fa"
with open(filepath, mode = "r") as handle:
    for record in tqdm(SeqIO.parse(handle, 'fasta')):
        accession = record.id
        if accession in set(allaccessions):
            totalcounts += 1
            for f in range(0, len(files)):
                if accession in accessions[f]:
                    coveredacc.append(accession)
                    counts[f] += 1
                    outfiles[f].write("\n"+accession+"  ")
                    sequence = list(record.seq)[21562:25384]
                    sequence = "".join(sequence)
                    outfiles[f].write(sequence+"\n")

print(counts)
print(totalcounts)

for x in allaccessions:
    if x not in coveredacc:
        print(x)



treefiles = os.listdir("trees")
#treefiles.remove(".DS_Store")

for f in treefiles:
    tree = open("trees/"+f, "r").readlines()[0]
    outfile = open("trees/"+f, "w")
    outfile.write(" 1000 1\n")
    tree = re.sub("\\).\\....:", "):", tree)
    tree = re.sub("\\).\\....;", ");", tree)
    outfile.write(tree)

print(outfilenames)
print(treefiles)
for i in range(0, len(files)):
    infile = open("basemltemplate.ctl", "r").readlines()
    outfile = open("baseml.ctl", "w")
    for line in infile:
        outfile.write(line.replace("SEQUENCES", outfilenames[i]).replace("TREEFILE", treefiles[i]))
    outfile.close()
    try:
        subprocess.call("paml/bin/baseml", shell = True)
        subprocess.call("Rscript harvestrates.R", shell = True)
    except:
        continue


