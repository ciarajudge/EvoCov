import os
import sys
import subprocess
from tqdm import tqdm
from operator import add
import numpy as np

form = sys.argv[1]
infile = open(sys.argv[2], "r").readlines()
meta = sys.argv[3]
targets = sys.argv[5]
if meta == "location":
    m = 3
elif meta == "variant":
    m = 2
else:
    m = 1

outfile = open(sys.argv[4], "w")

count = 0
for line in tqdm(infile):
    if list(line)[0] == ">":
        posns = []
        factor = line.split("|")[m]
        if factor == targets:
            count += 1
            counting = True
            outfile.write(line)
        else:
            counting = False
    else:
        if counting == True:
            outfile.write(line)
        else:
            continue


print(count)
