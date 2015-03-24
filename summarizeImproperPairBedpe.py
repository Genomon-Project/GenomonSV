#! /usr/local/bin/python

"""
    script for merging and gummarizing improper pair information 
"""

import sys

inputFile = sys.argv[1]

checkMarginSize = 1200

hIN = open(inputFile, 'r')

mergedBedpe = {}
for line in hIN:

    F = line.rstrip('\n').split('\t')

    # if F[6] == "HWI-ST1021:119:C14DVACXX:3:2307:12722:83115":
    #     print F[6]

    match = 0
    delList = []
    for key in sorted(mergedBedpe):

        tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2 = key.split('\t')
        tids, tmqs, talns = mergedBedpe[key].split('\t')

        # if F[6] == "HWI-ST1021:119:C14DVACXX:2:1203:19326:189526":
        #     print tids

        if F[0] != tchr1 or int(F[1]) > int(tend1) + checkMarginSize:
            talns_a = talns.split(';')
            talns_a_uniq = list(set(talns_a))

            if len(talns_a_uniq) >= 1:
                
                print '\t'.join([tchr1, tstart1, tend1, tchr2, tstart2, tend2, \
                                 tids, tmqs, tdir1, tdir2, talns])
                delList.append(key)
                continue

        else:
            
            if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2:
                if int(F[2]) > int(tstart1) and int(F[1]) <= int(tend1) and int(F[5]) > int(tstart2) and int(F[4]) <= int(tend2):

                    match = 1
                    newStart1 = str(max(int(tstart1), int(F[1])))  
                    newEnd1 = str(min(int(tend1), int(F[2])))
                    newStart2 = str(max(int(tstart2), int(F[4])))
                    newEnd2 = str(min(int(tend2), int(F[5])))

                    newKey = '\t'.join([tchr1, newStart1, newEnd1, tchr2, newStart2, newEnd2, tdir1, tdir2])
                    newIds = tids + ';' + F[6]
                    newMqs = tmqs + ';' + F[7]
                    newAlns = talns + ';' + F[10]
    
                    if newKey != key:    
                        delList.append(key)

                    mergedBedpe[newKey] = newIds + '\t' + newMqs + '\t' + newAlns
                    break

    for item in delList:
        del mergedBedpe[item]

    if match == 0:
        newKey = '\t'.join([F[0], F[1], F[2], F[3], F[4], F[5], F[8], F[9]])
        mergedBedpe[newKey] = F[6] + '\t' + F[7] + '\t' + F[10]


hIN.close()


for key in sorted(mergedBedpe):

    tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2 = key.split('\t')
    tids, tmqs, talns = mergedBedpe[key].split('\t')

    talns_a = talns.split(';')
    talns_a_uniq = list(set(talns_a))

    if len(talns_a_uniq) >= 1:

        print '\t'.join([tchr1, tstart1, tend1, tchr2, tstart2, tend2, \
                         tids, tmqs, tdir1, tdir2, talns])


