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

for line in hIN:
    F = line.rstrip('\n').split('\t')

    controlFlag = 0
    records = control_tb.query(F[0], int(F[1]), int(F[1]) + 1)
    for record in records:
        if "\t".join(F[0:6]) == "\t".join(record[0:6]):
            controlSamples = record[6].split(',')
            controlNums = record[7].split(',')
            for i in range(0, len(controlSamples)):
                if controlSamples[i] != matchedNormal is not None and int(controlNums[i]) >= int(supportReadThres):
                # if controlSamples[i] != matchedNormal and int(controlNums[i]) >= int(supportReadThres):
                    controlFlag = 1
            
    if controlFlag == 0:
        print "\t".join(F)

hIN.close()

