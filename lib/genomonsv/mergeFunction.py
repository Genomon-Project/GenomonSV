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

        

def organizeControl(inputFilePath, outputFilePath, check_margin_size):

    """
        script for organizing control junction information
    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    mergedBedpeInfo = {}
    mergedreadNums = {}
    num = 1
    for line in hIN:

        F = line.rstrip('\n').split('\t')

        match = 0
        delList = []
        for key in sorted(mergedBedpeInfo):

            tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2, inseqSize = key.split('\t')
            tsamples, treadNums = mergedBedpeInfo[key].split('\t')

            if tend1 == "18606952": print key

            # the investigated key is sufficiently far from the current line in the input file and no additional line to merge is expected. therefore flush the key and information
            if F[0] != tchr1 or int(F[2]) > int(tend1) + check_margin_size:

                # treatment for flush!!
                print >> hOUT, '\t'.join([tchr1, tstart1, tend1, tchr2, tstart2, tend2, "controlJunction_" + str(num), inseqSize, tdir1, tdir2, tsamples, treadNums])
                num = num + 1

                delList.append(key)
                continue


            else:

                # check whether the investigated key and the current line should be merged or not 
                if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2:

                    flag = 0
                    # detailed check on the junction position considering inserted sequences
                    if F[8] == "+":
                        expectedDiffSize = (int(F[2]) - int(tend1)) + (int(F[7]) - int(inseqSize))
                        if (F[9] == "+" and int(F[5]) == int(tend2) - int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(tend2) + int(expectedDiffSize)):
                            flag = 1
                    else:
                        expectedDiffSize = (int(F[2]) - int(tend1)) + (int(inseqSize) - int(F[7]))
                        if (F[9] == "+" and int(F[5]) == int(tend2) + int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(tend2) - int(expectedDiffSize)):
                            flag = 1

                    # if the junction position and direciton match
                    if flag == 1:

                        match = 1
                        newKey = '\t'.join([F[0], F[1], F[2], F[3], F[4], F[5], F[8], F[9], F[7]])
                        newSamples = tsamples + ';' + F[10]
                        newReadNums = treadNums + ';' + F[11]
                    
                        if F[7] < inseqSize:
                            mergedBedpeInfo[newKey] = '\t'.join([newSamples, newReadNums])
                            delList.append(key)
                        else:
                            mergedBedpeInfo[key] = '\t'.join([newSamples, newReadNums])

        for item in delList:
            del mergedBedpeInfo[item]

        # if the current line in the input file does not match any of the pooled keys
        if match == 0:
            newKey = '\t'.join([F[0], F[1], F[2], F[3], F[4], F[5], F[8], F[9], F[7]])
            mergedBedpeInfo[newKey] = F[10] + '\t' + F[11]



    hIN.close()


    for key in sorted(mergedBedpeInfo):

        tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2, inseqSize = key.split('\t')
        tsamples, treadNums = mergedBedpeInfo[key].split('\t')

        # treatment for flush!!
        print >> hOUT, '\t'.join([tchr1, tstart1, tend1, tchr2, tstart2, tend2, "controlJunction_" + str(num), inseqSize, tdir1, tdir2, tsamples, treadNums])
        num = num + 1


    hOUT.close()


 
