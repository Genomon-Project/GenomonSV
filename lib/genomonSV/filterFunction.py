#!/usr/bin/env python

"""
     functions for filtering candidates of structural variations
"""

import sys, gzip, pysam

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




def filterNonMatchControl(inputFilePath, outputFilePath, controlFile, matchedNormal, Params):

    """
        script for removing candidate in which 
        non-matched normals have the junction reads
    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    tabixfile = pysam.TabixFile(controlFile)

    controlPanel_num_thres = Params["controlPanel_num_thres"]
    controlPanel_check_margin = Params["controlPanel_check_margin"]

    
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        inseqSize = (0 if F[7] == "---" else len(F[7]))

        controlFlag = 0

        ####################
        # get the records for control junction data for the current position
        tabixErrorFlag = 0
        try:
            records = tabixfile.fetch(F[0], int(F[1]) - controlPanel_check_margin, int(F[2]) + controlPanel_check_margin)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1
        ####################


        ####################
        # for each record in control junction extracted, check the consistency with the current junction
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')

                if F[0] == record[0] and F[3] == record[3] and F[8] == record[8] and F[9] == record[9]:

                    flag = 0
                    # detailed check on the junction position considering inserted sequences
                    if F[8] == "+":
                        expectedDiffSize = (int(F[2]) - int(record[2])) + (inseqSize - int(record[7]))
                        if (F[9] == "+" and int(F[5]) == int(record[5]) - int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(record[5]) + int(expectedDiffSize)):
                            flag = 1
                    else:
                        expectedDiffSize = (int(F[2]) - int(record[2])) + (int(record[7]) - inseqSize)
                        if (F[9] == "+" and int(F[5]) == int(record[5]) + int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(record[5]) - int(expectedDiffSize)):
                            flag = 1

                    # if position relationship including inserted sequences matches
                    if flag == 1:
                        controlSamples = record[10].split(';')
                        controlNums = record[11].split(';')
                        for i in range(0, len(controlSamples)):
                            if controlSamples[i] != matchedNormal is not None and int(controlNums[i]) >= int(controlPanel_num_thres):
                            # if controlSamples[i] != matchedNormal and int(controlNums[i]) >= int(supportReadThres):
                                controlFlag = 1
        ####################
                    
        if controlFlag == 0:
            print >> hOUT, "\t".join(F)

    hIN.close()
    hOUT.close()
    tabixfile.close()



def addImproperInfo(inputFilePath, outputFilePath, improperFilePath):

    juncFile = sys.argv[1]
    improperFile = sys.argv[2]

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')
    tabixfile = pysam.TabixFile(improperFilePath)

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        improper_readNames = "---"
        improper_MQs = "---"
        improper_coveredRegions = "---"

        # get the records for control junction data for the current position
        tabixErrorFlag = 0
        try:
            records = tabixfile.fetch(F[0], int(F[1]) - -1, int(F[2]) + 1)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1


        if tabixErrorFlag == 0:
            for record_line in records:
                FF = record_line.split('\t')

                if F[0] == FF[0] and F[3] == FF[3] and int(F[1]) >= int(FF[1]) and int(F[2]) <= int(FF[2]) and int(F[4]) >= int(FF[4]) and int(F[5]) <= int(FF[5]) and F[8] == FF[8] and F[9] == FF[9]:
                    improper_readNames = FF[6]
                    improper_MQs = FF[7]
                    improper_coveredRegions = FF[10]
    

        print >> hOUT, "\t".join(F) + "\t" + improper_readNames + "\t" + improper_MQs + "\t" + improper_coveredRegions

    hIN.close()
    hOUT.close()
    tabixfile.close()


