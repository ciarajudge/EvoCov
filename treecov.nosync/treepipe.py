import os
import sys
import re
import subprocess
from Bio import SeqIO
from tqdm import tqdm
from operator import add

success = 0

try:
    for rep in range(0, 1):
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
        numseqs = len(accs)

        outfiles = []
        outfilenames = []
        for f in range(0, len(files)):
            outfiles.append(open("sequences/seqs"+str(f+1)+".nuc", "w"))
            outfilenames.append("sequences/seqs"+str(f+1)+".nuc")
            outfiles[f].write(" "+str(numseqs)+" 3822\n\n")

        counts = [0]*len(files)
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
            outfile.write(" "+str(numseqs)+" 1\n")
            tree = re.sub("\\).\\....:", "):", tree)
            tree = re.sub("\\).\\....;", ");", tree)
            outfile.write(tree)

        print(outfilenames)
        print(treefiles)
        for i in range(0, len(files)):
            if counts[i] < numseqs:
                continue
            success += 1
            infile = open("basemltemplate.txt", "r").readlines()
            outfile = open("baseml"+str(i)+".ctl", "w")
            for line in infile:
                outfile.write(line.replace("SEQUENCES", "sequences/seqs"+str(i+1)+".nuc").replace("TREEFILE", "sampledtree"+str(i+1)+".nhx"))
            outfile.close()
            '''
            try:
                x = input("BASEML TIME")
                if x == "success":
                    subprocess.call("Rscript harvestrates.R", shell = True)
            except:
                continue
            '''

    print(success)
    
except:
    subprocess.call("curl http://textbelt.com/text -d number=353877910680 -d message=\"Somethings Gone Wrong\" -d key=e237256cdcb7af72df888a7558f92c0e97b0fb55OVbGpYEsvYujZXDCVi0Rtvom6 ", shell = "True")
    
