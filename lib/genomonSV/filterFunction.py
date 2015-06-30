#!/usr/bin/env python

"""
     functions for filtering candidates of structural variations
"""

import sys, gzip, pysam, numpy
import coveredRegions

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



def filterMergedJunc(inputFilePath, outputFilePath, Params):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    min_support_num = Params["min_support_num"]
    min_mapping_qual = Params["min_mapping_qual"]
    min_cover_size = Params["min_cover_size"]


    for line in hIN:
        F = line.rstrip('\n').split('\t')

        ##########
        MQs = F[10].split(';')
        coveredRegion_primary = F[11].split(';')
        coveredRegion_pair = F[13].split(';')
        pairMQs = map(lambda x: int(x), F[12].split(';'))
        juncPairPos = F[14].split(';')
        juncReadTypes = F[15].split(';')
        improperCoveredRegion = ([] if F[19] == "---" else F[19].split(';'))

        # enumerate support read number
        juncSupport = len(MQs)
        improperMQs = ([] if F[18] == "---" else F[18].split(';'))
        improperSupport = len(improperMQs)
        # skip if the number of suppor read is below the minSupReadNum 
        if juncSupport + improperSupport < min_support_num:
            continue


        ####################
        # check for the mapping quality
        MQs1 = [];
        MQs2 = [];
        for i in range(0, len(MQs)):
            if juncReadTypes[i] == "1":
                MQs1.append(int(MQs[i]))
            else:
                MQs2.append(int(MQs[i]))

        # improper pair
        IMQs1 = [];
        IMQs2 = [];
        for i in range(0, len(improperMQs)):
            tempMQ = improperMQs[i].split(',')
            IMQs1.append(int(tempMQ[0]))
            IMQs2.append(int(tempMQ[1]))

     
        mapFlag = 0
        if len(MQs1) > 0 and numpy.median(MQs1) >= min_mapping_qual:
            mapFlag = 1

        if len(MQs2) > 0 and numpy.median(MQs2) >= min_mapping_qual:
            mapFlag = 1     

        if len(IMQs1) > 0 and numpy.median(IMQs1) >= min_mapping_qual:
            mapFlag = 1

        if len(IMQs2) > 0 and numpy.median(IMQs2) >= min_mapping_qual:
            mapFlag = 1

        if mapFlag == 0:
            continue


        #################### 
        # check for the covered region
        region1 = coveredRegions.coveredRegions()
        region2 = coveredRegions.coveredRegions()

        for i in range(0, len(coveredRegion_primary)):
            regPair = coveredRegion_primary[i].split(',')
            if juncReadTypes[i] == "1":
                region1.addMerge(regPair[0])
                region2.addMerge(regPair[1])
            else:
                region1.addMerge(regPair[1])
                region2.addMerge(regPair[0])


        for i in range(0, len(coveredRegion_pair)):
            if (juncReadTypes[i] == "1" and juncPairPos[i] == "1") or (juncReadTypes[i] == "2" and juncPairPos[i] == "2"): 
                region1.addMerge(coveredRegion_pair[i])
            elif (juncReadTypes[i] == "1" and juncPairPos[i] == "2") or (juncReadTypes[i] == "2" and juncPairPos[i] == "1"):
                region2.addMerge(coveredRegion_pair[i])


        for i in range(0, len(improperCoveredRegion)):
            regPair = improperCoveredRegion[i].split(',')
            region1.addMerge(regPair[0])
            region2.addMerge(regPair[1])


        region1.reduceMerge()
        region2.reduceMerge()

        if region1.regionSize() < min_cover_size or region2.regionSize() < min_cover_size:
            continue

        # every condition is satisfied
        print >> hOUT, '\t'.join(F)


    hIN.close()
    hOUT.close()

