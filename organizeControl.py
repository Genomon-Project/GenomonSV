#! /usr/local/bin/python

"""
    script for organizing control junction information
"""

import sys

inputFile = sys.argv[1]

hIN = open(inputFile, 'r')

tempKey = ""
tempSamples = ""
tempNums = ""
for line in hIN:
    F = line.rstrip('\n').split('\t')
    key = '\t'.join(F[0:6]) + '\t' + F[8] + '\t' + F[9]

    if key != tempKey:

        if tempKey != "":
            FF = key.split('\t')
            print '\t'.join(FF[0:6]) + '\t' + tempSamples + '\t' + tempNums + '\t' + FF[6] + '\t' + FF[7]

        tempKey = key
        tempSamples = F[6]
        tempNums = F[7]

    else:
        tempSamples = tempSamples + ',' + F[6]
        tempNums = tempNums + ',' + F[7]

hIN.close()


if tempKey != "":
    FF = key.split('\t')
    print '\t'.join(FF[0:6]) + tempSamples + '\t' + tempNums + '\t' + FF[6] + FF[7]

