#! /usr/local/bin/python

"""
    script for filtering by the length of SV size and support read junctions
"""

import sys, gzip

inputFile = sys.argv[1]
minJuncNum = sys.argv[2]
minSize = sys.argv[3]

hIN = gzip.open(inputFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')

    ##########
    # for now only consider long ranged SV??
    if F[0] == F[3] and abs(int(F[2]) - int(F[5])) < int(minSize): continue 
    ##########


    MQs = F[7].split(';')

    # enumerate support read number
    juncSupport = int(len(MQs))

    # skip if the number of suppor read is below the minSupReadNum 
    if juncSupport < int(minJuncNum): continue

    print '\t'.join(F)

hIN.close()
