#! /usr/local/bin/perl

import sys
import tabix

juncFile = sys.argv[1]
improperFile = sys.argv[2]

hIN = open(juncFile, 'r')
improper_tb = tabix.open(improperFile)

for line in hIN:
    F = line.rstrip('\n').split('\t')

    improper_readNames = "---"
    improper_MQs = "---"
    improper_coveredRegions = "---"

    tabixErrorFlag = 0
    try:
        records = improper_tb.query(F[0], int(F[1]), int(F[1]) + 1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for FF in records:

            if F[0] == FF[0] and F[3] == FF[3] and int(F[1]) >= int(FF[1]) and int(F[2]) <= int(FF[2]) and int(F[4]) >= int(FF[4]) and int(F[5]) <= int(FF[5]) and F[8] == FF[8] and F[9] == FF[9]:
                improper_readNames = FF[6]
                improper_MQs = FF[7]
                improper_coveredRegions = FF[10]
    

    print "\t".join(F) + "\t" + improper_readNames + "\t" + improper_MQs + "\t" + improper_coveredRegions

hIN.close()

