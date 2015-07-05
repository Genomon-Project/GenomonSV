#! /usr/bin/env python

import gzip


def simplifyJunc(inputFilePath, outputFilePath, label):

    """
        function for creating control junction information for filtering
    """

    hIN = gzip.open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'a')
    
    num = 1
    for line in hIN:

        F = line.strip('\n').split('\t')
        MQs = F[10].split(';')
        inseqLen = 0

        if F[7] != "---": inseqLen = len(F[7])

        print >> hOUT, '\t'.join(F[0:6]) + '\t' + "junction_" + str(num)  + '\t' + str(inseqLen) + '\t' + F[8] + '\t' + F[9] + '\t' + label + '\t' + str(len(MQs))
        num = num + 1

    hIN.close()
    hOUT.close()

         
