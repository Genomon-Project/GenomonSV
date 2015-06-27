#!/usr/bin/env python

"""
     functions for filtering candidates of structural variations
"""

import gzip


def filterJuncNumAndSize(inputFilePath, outputFilePath, Params):
     
    """
        script for filtering by the length of SV size and support read junctions
    """

    junc_num_thres = Params["junc_num_thres"]
    SV_size_thres = Params["SV_size_thres"]

    hIN = gzip.open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        ##########
        # for now only consider long ranged SV??
        svLen = abs(int(F[2]) - int(F[5]))
        if F[7] != "---": svLen = svLen + len(F[7])
        if F[0] == F[3] and svLen < int(SV_size_thres): continue
        ########## 
    
        ##########
        # enumerate support read number
        IDs = F[6].split(';')
        juncSupport = int(len(IDs))

        # skip if the number of suppor read is below the minSupReadNum 
        if juncSupport < int(junc_num_thres): continue
        ##########

        print >> hOUT, '\t'.join(F)


    hIN.close()
    hOUT.close()


