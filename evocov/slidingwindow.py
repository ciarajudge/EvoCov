import sys
import os
import numpy as np
import csv
import math

outfile = open("slidingwindowcandidates.txt", "w")

for n in range(0, 536):
    for x in range(0,7):
        outfile.write(str(n+x)+",")
    outfile.write("_3\n")
