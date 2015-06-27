#! /usr/local/bin/python

"""
    script for creating control junction information for filtering
"""

import sys, gzip

inputFile = sys.argv[1]
label = sys.argv[2]

num = 1
hIN = gzip.open(inputFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')

    MQs = F[10].split(';')
    inseqLen = 0
    if F[7] != "---": inseqLen = len(F[7])

    print '\t'.join(F[0:6]) + '\t' + "junction_" + str(num)  + '\t' + str(inseqLen) + '\t' + F[8] + '\t' + F[9] + '\t' + label + '\t' + str(len(MQs))
    num = num + 1

hIN.close() 
    
