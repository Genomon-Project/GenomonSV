#! /usr/local/bin/python

"""
    script for merging and summarizing junction read pairs
"""

import sys

inputFile = sys.argv[1]

checkMarginSize = 1000


hIN = open(inputFile, 'r')

mergedBedpe = {}
for line in hIN:

    F = line.rstrip('\n').split('\t')

    match = 0
    delList = []
    for key in sorted(mergedBedpe):

        tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2, inseqSize = key.split('\t')
        tids, tmqs1, talns1, tinseqs, tmqs2, talns2, tpinds, tcinds = mergedBedpe[key].split('\t') 

        if F[0] != tchr1 or int(F[1]) > int(tend1) + checkMarginSize:

            talns_a1 = talns1.split(';')
            talns_a1_uniq = list(set(talns_a1))

            talns_a2 = talns2.split(';')
            talns_a2_uniq = list(set(talns_a2))

            if len(talns_a2_uniq) >= 1:

                print '\t'.join([tchr1, tstart1, tend1, tchr2, tstart2, tend2, \
                                 tids, tmqs1, tdir1, tdir2, talns1, tinseqs, \
                                 tmqs2, talns2, tpinds, tcinds]) 
                delList.append(key)
                continue

        else:
 
            if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2:

                flag = 0
                if F[8] == "+":
                    expectedDiffSize = (int(F[2]) - int(tend1)) + (len(F[11]) - int(inseqSize))
                    if (F[9] == "+" and int(F[5]) == int(tend2) - int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(tend2) + int(expectedDiffSize)):
                        flag = 1
                else:
                    expectedDiffSize = (int(F[2]) - int(tend1)) + (int(inseqSize) - len(F[11]))
                    if (F[9] == "+" and int(F[5]) == int(tend2) + int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(tend2) - int(expectedDiffSize)):
                        flag = 1

                if flag == 1:
            
                    match = 1
                    newIds = tids + ';' + F[6]
                    newMqs1 = tmqs1 + ';' + F[7]
                    newAlns1 = talns1 + ';' + F[10]
                    newInseqs = tinseqs + ';' + F[11]
                    newPinds = tpinds + ';' + F[14]
                    newCinds = tcinds + ';' + F[15]
                    newMqs2 = tmqs2 + ';' + F[12]
                    newAlns2 = talns2 + ';' + F[13]


                    mergedBedpe[key] = '\t'.join([newIds, newMqs1, newAlns1, newInseqs, newMqs2, newAlns2, newPinds, newCinds])

    for item in delList:
        del mergedBedpe[item]

    if match == 0:
        newKey = '\t'.join([F[0], F[1], F[2], F[3], F[4], F[5], F[8], F[9], str(len(F[11]))])
        mergedBedpe[newKey] = F[6] + '\t' + F[7] + '\t' + '\t'.join(F[10:16])

hIN.close()

