#! /usr/local/bin/perl

"""
    script for removing candidate in which 
    non-matched normals have the junction reads
"""

import sys
import tabix
import re

inputFile = sys.argv[1]
controlFile = sys.argv[2]
matchedNormal = sys.argv[3]
supportReadThres = sys.argv[4]

hIN = open(inputFile, 'r')
control_tb = tabix.open(controlFile)
margin = 40 # the positions of junction slip very large (sometimes ove 20bp). I think this value should be at least 30bp, but may be more) 

for line in hIN:
    F = line.rstrip('\n').split('\t')

    inseqSize = (0 if F[7] == "---" else len(F[7]))

    controlFlag = 0
    records = control_tb.query(F[0], int(F[1]) - margin, int(F[2]) + margin)
    for record in records:

        if F[0] == record[0] and F[3] == record[3] and F[8] == record[8] and F[9] == record[9]:

            flag = 0
            # detailed check on the junction position considering inserted sequences
            if F[8] == "+":
                expectedDiffSize = (int(F[2]) - int(record[2])) + (inseqSize - int(record[7]))
                if (F[9] == "+" and int(F[5]) == int(record[5]) - int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(record[5]) + int(expectedDiffSize)):
                    flag = 1
            else:
                expectedDiffSize = (int(F[2]) - int(record[2])) + (int(record[7]) - inseqSize)
                if (F[9] == "+" and int(F[5]) == int(record[5]) + int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(record[5]) - int(expectedDiffSize)):
                    flag = 1

            # if position relationship including inserted sequences matches
            if flag == 1:
                controlSamples = record[10].split(';')
                controlNums = record[11].split(';')
                for i in range(0, len(controlSamples)):
                    if controlSamples[i] != matchedNormal is not None and int(controlNums[i]) >= int(supportReadThres):
                    # if controlSamples[i] != matchedNormal and int(controlNums[i]) >= int(supportReadThres):
                        controlFlag = 1
            
    if controlFlag == 0:
        print "\t".join(F)

hIN.close()

