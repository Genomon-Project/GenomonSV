#!/usr/bin/env python

"""
     functions for filtering candidates of structural variations
"""

import sys, gzip, subprocess, pysam, numpy, math, os, re
import coveredRegions
import realignmentFunction
import utils
from scipy import stats

def filterJuncNumAndSize(inputFilePath, outputFilePath, junc_num_thres, sv_size_thres, inversion_size_thres):
     
    """
        script for filtering by the length of SV size and support read junctions
    """

    hIN = gzip.open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        ##########
        # for now only consider long ranged SV??
        svLen = abs(int(F[2]) - int(F[5]))
        if F[7] != "---": svLen = svLen + len(F[7])
        if F[0] == F[3]:
            # for inversion, use other threshould 
            if F[8] == F[9] and svLen < int(inversion_size_thres):
                continue
            if svLen < int(sv_size_thres): continue
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




def filterNonMatchControl(inputFilePath, outputFilePath, controlFile, matchedNormal, controlPanel_num_thres, controlPanel_check_margin):

    """
        script for removing candidate in which 
        non-matched normals have the junction reads
    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    use_control = True if controlFile != "" else False
    if use_control == True: tabixfile = pysam.TabixFile(controlFile)

    tabixErrorMsg = "" 
    for line in hIN:
        F = line.rstrip('\n').split('\t')


        controlFlag = 0
        max_control_sample = "---"
        max_control_num = 0

        if use_control == True:

            inseqSize = (0 if F[7] == "---" else len(F[7]))

            ####################
            # get the records for control junction data for the current position
            tabixErrorFlag = 0
            try:
                records = tabixfile.fetch(F[0], int(F[1]) - controlPanel_check_margin, int(F[2]) + controlPanel_check_margin)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorMsg = str(inst.args)
                tabixErrorFlag = 1
            ####################


            ####################
            # for each record in control junction extracted, check the consistency with the current junction
            # max_control_sample = "---"
            # max_control_num = 0 
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
                                if controlSamples[i] == matchedNormal: continue

                                if int(controlNums[i]) > max_control_num:
                                    max_control_sample = controlSamples[i]
                                    max_control_num = int(controlNums[i])

                                if int(controlNums[i]) >= int(controlPanel_num_thres):
                                    controlFlag = 1

                                """
                                # if controlSamples[i] != matchedNormal is not None and int(controlNums[i]) >= int(controlPanel_num_thres):
                                # if controlSamples[i] != matchedNormal and int(controlNums[i]) >= int(supportReadThres):
                                    controlFlag = 1
                                    if int(controlNums[i]) > max_control_num:
                                        max_control_sample = controlSamples[i]
                                        max_control_num = int(controlNums[i]) 
                                """
                            
            ####################
              
        if controlFlag == 0:
            print >> hOUT, "\t".join(F) + '\t' + max_control_sample + '\t' + str(max_control_num)


    if tabixErrorMsg != "":
        utils.warningMessage("One or more error occured in tabix file fetch, e.g.: " + tabixErrorMsg)

    hIN.close()
    hOUT.close()
    if use_control == True: tabixfile.close()



def addImproperInfo(inputFilePath, outputFilePath, improperFilePath):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')
    tabixfile = pysam.TabixFile(improperFilePath)

    tabixErrorMsg = ""
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        improper_readNames = "---"
        improper_MQs = "---"
        improper_coveredRegions = "---"

        # get the records for control junction data for the current position
        tabixErrorFlag = 0
        try:
            records = tabixfile.fetch(F[0], int(F[1]) - 1, int(F[2]) + 1)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1


        if tabixErrorFlag == 0:
            for record_line in records:
                FF = record_line.split('\t')

                if F[0] == FF[0] and F[3] == FF[3] and int(F[1]) >= int(FF[1]) and int(F[2]) <= int(FF[2]) and int(F[4]) >= int(FF[4]) and int(F[5]) <= int(FF[5]) and F[8] == FF[8] and F[9] == FF[9]:
                    improper_readNames = FF[6]
                    improper_MQs = FF[7]
                    improper_coveredRegions = FF[10]
    
        # last two columns are about pooled control information
        print >> hOUT, "\t".join(F[:-2]) + "\t" + improper_readNames + "\t" + improper_MQs + "\t" + improper_coveredRegions + '\t' + '\t'.join(F[-2:])

    if tabixErrorMsg != "":
        utils.warningMessage("One or more error occured in tabix file fetch, e.g.: " + tabixErrorMsg)


    hIN.close()
    hOUT.close()
    tabixfile.close()



def filterMergedJunc(inputFilePath, outputFilePath, min_support_num, min_mapping_qual, min_cover_size):

    """    
    script for filtering the inferred break points for structural variation

    as conditions for filtering
    min_support_num: the number of support reads
    min_mapping_qual: (e.g., in either type 1 or 2 junction reads, median of mapping quality is at least 40 ??)
    min_cover_size: the size of alignment covered region (e.g., around the both break points, at least 80 bp sequences have to be covered)
    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')


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
        junc_ids = [re.sub(r'/\d$', '', x) for x in F[6].split(';')]
        # improper_ids = F[17].split(';')
        improper_ids = [] if F[17] == "---" else F[17].split(';')
        if len(list(set(junc_ids + improper_ids))) < min_support_num:
            continue

        juncSupport = len(MQs)
        improperMQs = ([] if F[18] == "---" else F[18].split(';'))
        """
        improperSupport = len(improperMQs)
        # skip if the number of suppor read is below the minSupReadNum 
        if juncSupport + improperSupport < min_support_num:
            continue
        """
        

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
        print >> hOUT, '\t'.join(F) + '\t' + str(region1.regionSize()) + '\t' + str(region2.regionSize())


    hIN.close()
    hOUT.close()



def removeClose(inputFilePath, outputFilePath, close_check_margin, close_check_thres):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    key2info = {}
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        delList = []
        skipFlag = 0
        for tkey in key2info:

            tchr1, tstart1, tend1, tchr2, tstart2, tend2, inseq, tdir1, tdir2 = tkey.split('\t')

            if F[0] != tchr1 or int(F[1]) > int(tend1) + close_check_margin:

                print >> hOUT, key2info[tkey]
                delList.append(tkey)
        
            else:

                if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2 and abs(int(F[2]) - int(tend1)) <= close_check_margin and abs(int(F[5]) - int(tend2)) <= close_check_margin:

                    infos = key2info[tkey].split('\t')
                    if len(F[6].split(';')) < len(infos[6].split(';')) and len(F[6].split(';')) < int(close_check_thres):
                        skipFlag = 1
                    elif len(F[6].split(';')) > len(infos[6].split(';')) and len(infos[6].split(';')) < int(close_check_thres):
                        delList.append(tkey)


        for tkey in delList:
            del key2info[tkey]

        if skipFlag == 0:
            key2info['\t'.join(F[0:6] + F[7:10])] = '\t'.join(F)

    hIN.close()


    for tkey in key2info:
        print >> hOUT, key2info[tkey]


    hOUT.close()



def validateByRealignment(inputFilePath, outputFilePath, tumorBamFilePath, normalBamFilePath, reference_genome, blat_option,
                          short_tandem_reapeat_thres, max_depth, search_length, search_margin, split_refernece_thres, validate_sequence_length):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')
    blat_cmds = ("blat " + blat_option).split(' ')

    num = 1
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        # if len(F) >= 24:
        chr1, pos1, dir1, chr2, pos2, dir2, juncSeq, max_control_sample, max_control_num, overhang1, overhang2 = F[0], F[2], F[8], F[3], F[5], F[9], F[7], F[20], F[21], F[22], F[23]
        # else:
        #     chr1, pos1, dir1, chr2, pos2, dir2, juncSeq = F[0], F[2], F[8], F[3], F[5], F[9], F[7]
        #     max_control_sample, max_control_num, overhang1, overhang2 = "---", "---", "---", "---"

        STDFlag = 0
        if chr1 == chr2 and int(pos2) - int(pos1) < short_tandem_reapeat_thres and dir1 == "-" and dir2 == "+": STDFlag = 1

        ####################
        # extract short reads from tumor sequence data around the candidate
        fRet = realignmentFunction.extractSVReadPairs(tumorBamFilePath, outputFilePath + ".tmp.tumor.fa", chr1, pos1, dir1, chr2, pos2, dir2, 
                                                      max_depth, search_length, search_margin)
        if fRet == 1: continue

        if normalBamFilePath != "":
            # extract short reads from matched-control sequence data around the candidate
            fRet = realignmentFunction.extractSVReadPairs(normalBamFilePath, outputFilePath + ".tmp.normal.fa", chr1, pos1, dir1, chr2, pos2, dir2,
                                                          max_depth, search_length, search_margin)
            if fRet == 1: continue
        ####################
        
        ####################
        # generate reference sequence and sequence containing the presumed variants
        if int(pos1) - validate_sequence_length < 0: continue
        if int(pos2) - validate_sequence_length < 0: continue
        realignmentFunction.getRefAltForSV(outputFilePath + ".tmp.refalt.fa", chr1, pos1, dir1, chr2, pos2, dir2, juncSeq,
                                           reference_genome, split_refernece_thres, validate_sequence_length)

        ####################
        # alignment tumor short reads to the reference and alternative sequences
        FNULL = open(os.devnull, 'w')
        subprocess.call(blat_cmds + [outputFilePath + ".tmp.refalt.fa", outputFilePath + ".tmp.tumor.fa", outputFilePath + ".tmp.tumor.psl"], 
                        stdout = FNULL, stderr = subprocess.STDOUT)
        
        ####################
        # alignment normal short reads to the reference and alternative sequences
        if normalBamFilePath != "":
            subprocess.call(blat_cmds + [outputFilePath + ".tmp.refalt.fa", outputFilePath + ".tmp.normal.fa", outputFilePath + ".tmp.normal.psl"],
                            stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        ####################
        # summarize alignment results
        tumorRef, tumorAlt = realignmentFunction.summarizeRefAlt(outputFilePath + ".tmp.tumor.psl", STDFlag)

        normalRef, normalAlt = "---", "---"
        if normalBamFilePath != "":
            normalRef, normalAlt = realignmentFunction.summarizeRefAlt(outputFilePath + ".tmp.normal.psl", STDFlag)

        # fisher test
        lpvalue = "---"
        if normalBamFilePath != "":
            oddsratio, pvalue = stats.fisher_exact([[tumorRef, tumorAlt], [normalRef, normalAlt]], 'less')
            if pvalue < 1e-100: pvalue = 1e-100
            lpvalue = (- math.log(pvalue, 10) if pvalue < 1 else 0)
            lpvalue = str(round(lpvalue, 4)) 

        print >> hOUT, '\t'.join([chr1, pos1, dir1, chr2, pos2, dir2, juncSeq, \
                                  str(tumorRef), str(tumorAlt), str(normalRef), str(normalAlt), str(lpvalue), \
                                  max_control_sample, max_control_num, overhang1, overhang2])

        if num % 100 == 0:        
            print >> sys.stderr, "finished checking " + str(num) + " candidates"
        num = num + 1

    
    if num > 1:
        subprocess.call(["rm", outputFilePath + ".tmp.tumor.fa"])
        subprocess.call(["rm", outputFilePath + ".tmp.refalt.fa"])
        subprocess.call(["rm", outputFilePath + ".tmp.tumor.psl"])

        if normalBamFilePath != "":
            subprocess.call(["rm", outputFilePath + ".tmp.normal.fa"])
            subprocess.call(["rm", outputFilePath + ".tmp.normal.psl"])

    hIN.close()
    hOUT.close()



def filterNumAFFis(inputFilePath, outputFilePath, normalBamFilePath, 
                   min_tumor_read_pair, min_tumor_alleleFreq, max_control_read_pair, max_control_allele_freq, max_fisher_pvalue):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')
  
 
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        tumorAF = 0 
        if float(F[7]) + float(F[8]) > 0: tumorAF = float(F[8]) / (float(F[7]) + float(F[8]))     
        tumorAF = str(round(tumorAF, 4))

        normalAF = "---"
        if normalBamFilePath != "":
            normalAF = 0
            if float(F[9]) + float(F[10]) > 0: normalAF = float(F[10]) / (float(F[9]) + float(F[10]))
            normalAF = str(round(normalAF, 4))

        if int(F[8]) < int(min_tumor_read_pair): continue
        if float(tumorAF) < float(min_tumor_alleleFreq): continue

        if normalBamFilePath != "":
            if int(F[10]) > int(max_control_read_pair): continue
            if float(normalAF) > float(max_control_allele_freq): continue
            if 10**(-float(F[11])) > float(max_fisher_pvalue): continue

        print >> hOUT, '\t'.join(F[0:9]) + '\t' + tumorAF + '\t' + '\t'.join(F[9:11]) + '\t' + normalAF + '\t' + '\t'.join(F[11:])

    hIN.close()
    hOUT.close()


 
