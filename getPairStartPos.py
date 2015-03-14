#! /usr/local/bin/python

"""
    script for obtaining the position information about the pair read from the junction file 

"""

import sys, re

inputFile = sys.argv[1]
hIN = open(inputFile, 'r')

num = 1
reChrPos = re.compile('^(\w+):(\d+)')
for line in hIN:
    F = line.rstrip('\n').split('\t')
    ID = F[6]
    chr = ""
    pos = ""

    # obtain the information about the start site of the pair read 
    chrpos = reChrPos.search(F[12])
    if chrpos is not None:
        chr = chrpos.group(1)
        pos = chrpos.group(2)
    else:
        print '\t'.join(F)
        print "the 13th column did not match to (chr):(start) pattern"
        sys.exit()

    # change the pair read num
    if ID[-1:] == "1":
        ID = ID[:-1] + "2"
    else:
        ID = ID[:-1] + "1"

    print chr + '\t' + str(int(pos) - 1) + '\t' + pos + '\t' + ID + '\t' + str(num)

    num = num + 1

hIN.close()

