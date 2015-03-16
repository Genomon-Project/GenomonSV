#! /usr/local/bin/python

"""
    script for creating control junction information for filtering
"""

import sys

inputFile = sys.argv[1]
label = sys.argv[2]

hIN = open(inputFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')

    MQs = F[7].split(';')
    print '\t'.join(F[0:6]) + '\t' + label + '\t' + str(len(MQs)) + '\t' + F[8] + '\t' + F[9]

hIN.close() 
    
