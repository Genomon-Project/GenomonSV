#! /usr/local/bin/python

import sys, re

originalFile = sys.argv[1]
pairInfoFile = sys.argv[2]

hOriginalFile = open(originalFile, 'r')
hPairInfoFile = open(pairInfoFile, 'r')

line1 = hOriginalFile.readline().rstrip('\n')
line2 = hPairInfoFile.readline().rstrip('\n')
for line2 in hPairInfoFile:
    F2 = line2.rstrip('\n').split('\t')
    ID2 = F2[3]
    ID2 = re.sub(r'/\d$', '', ID2)

    line1 = hOriginalFile.readline()
    F1 = line1.rstrip('\n').split('\t')
    ID1 = F1[6]
    ID1 = re.sub(r'/\d$', '', ID1)
    
    while ID1 != ID2:
        print >> sys.stderr, "No pair information for %s" % ID1
        line1 = hOriginalFile.readline()
        F1 = line1.rstrip('\n').split('\t')
        ID1 = F1[6]
        ID1 = re.sub(r'/\d$', '', ID1)
    
    print '\t'.join(F1[0:12]) + '\t' + F2[5] + '\t' + F2[6] + '\t' + F1[13]

hOriginalFile.close()
hPairInfoFile.close()

