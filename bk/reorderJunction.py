#! /usr/local/bin/python

"""
    script just for changing the order of bedpe file

"""

import sys, re

inputFile = sys.argv[1]
hIN = open(inputFile, 'r')

for line in hIN:
    chr1, start1, end1, chr2, start2, end2, ID, mapQ, dir1, dir2, align, inseq, pairPos = line.rstrip('\n').split('\t')

    if (chr1 == "hs37d5" or chr2 == "hs37d5"): continue
    if (re.search(r'_s$', ID) is not None): continue
 
    if chr1 < chr2: 
        print '\t'.join([chr1, start1, end1, chr2, start2, end2, ID, mapQ, dir1, dir2, align, inseq, pairPos, "1"])
    elif chr1 > chr2:
        print '\t'.join([chr2, start2, end2, chr1, start1, end1, ID, mapQ, dir2, dir1, align, inseq, pairPos, "2"])
    else:
        if int(start1) <= int(start2):
            print '\t'.join([chr1, start1, end1, chr2, start2, end2, ID, mapQ, dir1, dir2, align, inseq, pairPos, "1"])
        else:
            print '\t'.join([chr2, start2, end2, chr1, start1, end1, ID, mapQ, dir2, dir1, align, inseq, pairPos, "2"])


hIN.close()

