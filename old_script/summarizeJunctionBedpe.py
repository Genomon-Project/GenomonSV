#! /usr/local/bin/python

"""
    script for merging and summarizing junction read pairs
"""

import sys
from Bio.Seq import Seq
from collections import Counter

inputFile = sys.argv[1]

checkMarginSize = 1000 # this value should be at least 1000. but very large value will lead to slow computation..


hIN = open(inputFile, 'r')

mergedBedpeInfo = {}
mergedJunction = {}
for line in hIN:

    F = line.rstrip('\n').split('\t')

    match = 0
    delList = []
    for key in sorted(mergedBedpeInfo):

        tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2, inseqSize = key.split('\t')
        tids, tinseqs, tmqs1, talns1, tmqs2, talns2, tpinds, tcinds = mergedBedpeInfo[key].split('\t') 

        # the investigated key is sufficiently far from the current line in the input file and no additional line to merge is expected. therefore flush the key and information
        if F[0] != tchr1 or int(F[1]) > int(tend1) + checkMarginSize:

            # obtain the most frequent junction
            junc_counter = Counter(mergedJunction[key].split(';'))
            best_junc = junc_counter.most_common(1)[0][0]
            btchr1, btend1, btdir1, btchr2, btend2, btdir2, btinseq = best_junc.split(',') 
            btstart1 = str(int(btend1) - 1)
            btstart2 = str(int(btend2) - 1)


            print '\t'.join([btchr1, btstart1, btend1, btchr2, btstart2, btend2, \
                             tids, btinseq, btdir1, btdir2, tmqs1, talns1, \
                             tmqs2, talns2, tpinds, tcinds]) + '\t' +  \
                  mergedJunction[key]

            # add to the deletion list (later the key will removed from the dictionaries)
            delList.append(key)
            continue

        else:

            # check whether the investigated key and the current line should be merged or not 
            if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2:

                flag = 0
                # detailed check on the junction position considering inserted sequences
                if F[8] == "+":
                    expectedDiffSize = (int(F[2]) - int(tend1)) + (len(F[7]) - int(inseqSize))
                    if (F[9] == "+" and int(F[5]) == int(tend2) - int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(tend2) + int(expectedDiffSize)):
                        flag = 1
                else:
                    expectedDiffSize = (int(F[2]) - int(tend1)) + (int(inseqSize) - len(F[7]))
                    if (F[9] == "+" and int(F[5]) == int(tend2) + int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(tend2) - int(expectedDiffSize)):
                        flag = 1

                # if the junction position and direciton match
                if flag == 1:
            
                    match = 1
                    newIds = tids + ';' + F[6]
                    newInseqs = tinseqs + ';' + F[7]
                    newMqs1 = tmqs1 + ';' + F[10]
                    newAlns1 = talns1 + ';' + F[11]
                    newMqs2 = tmqs2 + ';' + F[12]
                    newAlns2 = talns2 + ';' + F[13]
                    newPinds = tpinds + ';' + F[14]
                    newCinds = tcinds + ';' + F[15]


                    mergedBedpeInfo[key] = '\t'.join([newIds, newInseqs, newMqs1, newAlns1, newMqs2, newAlns2, newPinds, newCinds])

                    # check whether the inserted sequence should be reverse-complemented 
                    tinseq = F[7]
                    if F[7] != "---" and F[8] == F[9] and F[15] == "2":
                        tinseq = str(Seq(F[7]).reverse_complement())

                    mergedJunction[key] = mergedJunction[key] + ";" + ','.join([F[0], F[2], F[8], F[3], F[5], F[9], tinseq])

    for item in delList:
        del mergedBedpeInfo[item]
        del mergedJunction[item]

    # if the current line in the input file does not match any of the pooled keys
    if match == 0:
        newKey = '\t'.join([F[0], F[1], F[2], F[3], F[4], F[5], F[8], F[9], str(len(F[7]))])
        mergedBedpeInfo[newKey] = F[6] + '\t' + F[7] + '\t' + '\t'.join(F[10:16])

        # check whether the inserted sequence should be reverse-complemented
        tinseq = F[7] 
        if F[7] != "---" and F[8] == F[9] and F[15] == "2":
            tinseq = str(Seq(F[7]).reverse_complement())

        mergedJunction[newKey] = ','.join([F[0], F[2], F[8], F[3], F[5], F[9], tinseq])

hIN.close()


# last treatment
for key in sorted(mergedBedpeInfo):

    tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2, inseqSize = key.split('\t')
    tids, tinseqs, tmqs1, talns1, tmqs2, talns2, tpinds, tcinds = mergedBedpeInfo[key].split('\t')

    # obtain the most frequent junction
    junc_counter = Counter(mergedJunction[key].split(';'))
    best_junc = junc_counter.most_common(1)[0][0]
    btchr1, btend1, btdir1, btchr2, btend2, btdir2, btinseq = best_junc.split(',')
    btstart1 = str(int(btend1) - 1)
    btstart2 = str(int(btend2) - 1)

    print '\t'.join([btchr1, btstart1, btend1, btchr2, btstart2, btend2, \
                     tids, btinseq, btdir1, btdir2, tmqs1, talns1, \
                     tmqs2, talns2, tpinds, tcinds]) + '\t' +  \
          mergedJunction[key]

