#! /usr/local/bin/python

import sys

tumorFile = sys.argv[1]
normalFile = sys.argv[2]

key2num = {}
hIN = open(normalFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')
    key = '\t'.join(F[0:6]) + '\t' + F[8] + '\t' + F[9]

    # enumerate support read number
    MQs = F[7].split(';')
    juncSupport = len(MQs)
    improperSupport = (0 if F[16] == "---" else len(F[16].split(';')))

    key2num[key] = str(juncSupport) + '\t' + str(improperSupport)

hIN.close()

num = 1
hIN = open(tumorFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')
    key = '\t'.join(F[0:6]) + '\t' + F[8] + '\t' + F[9]

    # enumerate support read number
    MQs = F[7].split(';')
    juncSupport = len(MQs)
    improperSupport = (0 if F[16] == "---" else len(F[16].split(';')))

    normalInfo = (key2num[key] if key in key2num else "---" + '\t' + "---")

    print '\t'.join(F[0:6]) + '\t' + "genomonSV_" + str(num) + '\t' + "0" + '\t' + F[8] + '\t' + F[9] + '\t' + str(juncSupport) + '\t' + str(improperSupport) + '\t' + normalInfo

    num = num + 1

hIN.close()



