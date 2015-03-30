#! /usr/local/bin/python

import sys, gzip

inputFile = sys.argv[1]
hIN = gzip.open(inputFile, 'r')

for line in hIN:
    F = line.rstrip('\n').split('\t')

    chr = F[2]
    starts = F[9].split(',')
    ends = F[10].split(',')
    strand = F[3]
    exonNum = int(F[8])
    gene = F[1]
    symbol = F[12]


    for i in range(0, len(starts) - 1):
        key = chr + '\t' + starts[i] + '\t' + ends[i]
        if strand == "+":
            print key + '\t' + symbol + "(" + gene + ")" + "." + str(i) + '\t' + "0" + '\t' + "+"
        else:
            print key + '\t' + symbol + "(" + gene + ")" + "." + str(exonNum - i - 1) + '\t' + "0" + '\t' + "-"


hIN.close()


