#! /usr/local/bin/python

import sys, gzip

inputFile = sys.argv[1]
hIN = gzip.open(inputFile, 'r')

for line in hIN:
    F = line.rstrip('\n').split('\t')

    chr = F[2]
    start = F[4]
    end = F[7]
    strand = F[3]
    symbol = F[12]

    chr = chr.replace('chr', '')

    key = chr + '\t' + start + '\t' + end
    print key + '\t' + symbol + '\t' + "0" + '\t' + strand

hIN.close()


