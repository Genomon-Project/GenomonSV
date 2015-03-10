#! /usr/local/bin/perl

import sys
import tabix

inputFile = sys.argv[1]
controlFile = sys.argv[2]

inFile = open(inputFile, 'r')
control_tb = tabix.open(controlFile)

for line in inFile:
    F = line.rstrip('\n').split('\t')

    controlSample = "---"
    controlNum = "---"
    records = control_tb.query(F[0], int(F[1]), int(F[1]) + 1)
    for record in records:
        if "\t".join(F[0:6]) == "\t".join(record[0:6]):
            controlSample = record[6]
            controlNum = record[7]

    print "\t".join(F) + "\t" + controlSample + "\t" + controlNum


