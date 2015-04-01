#! /usr/local/bin/python

import sys, re, pysam, gzip, tabix, numpy, math
from scipy import stats


tumorInfoFile = sys.argv[1]
normalJunction = sys.argv[2]
normalImproper = sys.argv[3]
tumorBamPath = sys.argv[4]
normalBamPath = sys.argv[5]

##########
# function

# getting the "proper" read pairs covering the break point (considering some margin for judging "cover")
def getCoverReadsNum(bamfile, varIDs, chr_region, pos_region, searchLength, margin):

    readIDs = []    
    for read in bamfile.fetch(chr_region, int(pos_region), int(pos_region + searchLength)):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip if not proper pair
        if flags[1] == "0": continue

        # skip if either of the read pair is unmapped
        if flags[2] == "1" or flags[3] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        # skip if not reverse read
        if flags[4] == "0": continue

        # add some margins to define "cover"
        if read.pnext + margin <= pos_region <= read.aend - margin:
            if read.qname not in varIDs:
                readIDs.append(read.qname)

    return(readIDs)




hIN = open(tumorInfoFile, 'r')
normalJunction_tb = tabix.open(normalJunction)
normalImproper_tb = tabix.open(normalImproper)
tumorBam_ps = pysam.Samfile(tumorBamPath, 'rb')
normalBam_ps = pysam.Samfile(normalBamPath, 'rb')


SVNum = 1
margin = 40
for line in hIN:

    F = line.strip('\n').split('\t')
    inseqSize = (0 if F[7] == "---" else len(F[7]))

    # tumor junction read pair count
    tumorJunctionIDs_temp = F[6].split(';')
    tumorJunctionIDs = map(lambda x: re.sub(r'/\d$', '', x), tumorJunctionIDs_temp)    
 
    # tumor improper read pair count
    if F[17] != "---":
        tumorImproperIDs = F[17].split(';')
    else:
        tumorImproperIDs = []

    # tumor junction and improper read pair
    tumorJunctionImproperIDs = list(set(tumorJunctionIDs + tumorImproperIDs))


    # normal junction read pair count
    tabixErrorFlag = 0
    try:
        records = normalJunction_tb.query(F[0], int(F[1]) - margin, int(F[2]) + margin)
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    normalJunctionIDs = [];
    if tabixErrorFlag == 0:
        for record in records:

            if F[0] == record[0] and F[3] == record[3] and F[8] == record[8] and F[9] == record[9]:
    
                normalInseqSize = (0 if record[7] == "---" else len(record[7]))
                flag = 0
                # detailed check on the junction position considering inserted sequences
                if F[8] == "+":
                    expectedDiffSize = (int(F[2]) - int(record[2])) + (inseqSize - normalInseqSize)
                    if (F[9] == "+" and int(F[5]) == int(record[5]) - int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(record[5]) + int(expectedDiffSize)):
                        flag = 1
                else:
                    expectedDiffSize = (int(F[2]) - int(record[2])) + (normalInseqSize - inseqSize)
                    if (F[9] == "+" and int(F[5]) == int(record[5]) + int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(record[5]) - int(expectedDiffSize)):
                        flag = 1

               # if position relationship including inserted sequences matches
                if flag == 1:
                    normalJunctionIDs_temp1 = record[6].split(';')
                    normalJunctionIDs = map(lambda x: re.sub(r'/\d$', '', x), normalJunctionIDs_temp1)


    # normal improper read pair count
    tabixErrorFlag = 0
    try:
        records = normalImproper_tb.query(F[0], int(F[1]), int(F[1]) + 1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    normalImproperIDs = []
    if tabixErrorFlag == 0:
        for FF in records:
            if F[0] == FF[0] and F[3] == FF[3] and int(F[1]) >= int(FF[1]) and int(F[2]) <= int(FF[2]) and int(F[4]) >= int(FF[4]) and int(F[5]) <= int(FF[5]) and F[8] == FF[8] and F[9] == FF[9]:
                normalImproperIDs_temp = FF[6].split(';')
                if len(normalImproperIDs_temp) > len(normalImproperIDs):
                    normalImproperIDs = normalImproperIDs_temp

    # normal junction and improper read pair
    normalJunctionImproperIDs = list(set(normalJunctionIDs + normalImproperIDs))

    # print '\t'.join(F[0:6]) + '\t' + '\t'.join([str(len(tumorJunctionIDs)), str(len(tumorImproperIDs)), str(len(tumorJunctionImproperIDs)), \
    #                                             str(len(normalJunctionIDs)), str(len(normalImproperIDs)), str(len(normalJunctionImproperIDs))])

    # tumor surrounding read pair count
    # print F[0]
    tumorProperNum1 = getCoverReadsNum(tumorBam_ps, tumorJunctionImproperIDs, F[0], int(F[2]), 1000, 20)
    tumorProperNum2 = getCoverReadsNum(tumorBam_ps, tumorJunctionImproperIDs, F[3], int(F[5]), 1000, 20)
    tumorProperNum = list(set(tumorProperNum1 + tumorProperNum2))

    # normal surrounding read pair count
    normalProperNum1 = getCoverReadsNum(normalBam_ps, normalJunctionImproperIDs, F[0], int(F[2]), 1000, 20)
    normalProperNum2 = getCoverReadsNum(normalBam_ps, normalJunctionImproperIDs, F[3], int(F[5]), 1000, 20)
    normalProperNum = list(set(normalProperNum1 + normalProperNum2))

    # fisher test
    oddsratio, pvalue = stats.fisher_exact([[len(tumorJunctionImproperIDs), len(normalJunctionImproperIDs)], [len(tumorProperNum), len(normalProperNum)]])

    if pvalue < 1e-100:
        pvalue = 1e-100

    print '\t'.join(F[0:6]) + '\t' + "genomonSV_" + str(SVNum) + '\t' + F[7] + '\t' + F[8] + '\t' + F[9] + '\t' + \
          '\t'.join([str(len(tumorJunctionIDs)), str(len(tumorImproperIDs)), str(len(tumorJunctionImproperIDs)), str(len(tumorProperNum)), \
                     str(len(normalJunctionIDs)), str(len(normalImproperIDs)), str(len(normalJunctionImproperIDs)), str(len(normalProperNum))]) + '\t' + \
          str(- math.log(pvalue, 10))


    SVNum = SVNum + 1




