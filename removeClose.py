#! /usr/local/bin/python

import sys

inputFile = sys.argv[1]

hIN = open(inputFile, 'r')

key2info = {}
for line in hIN:
    F = line.rstrip('\n').split('\t')

    delList = []
    skipFlag = 0
    for tkey in key2info:

        tchr1, tstart1, tend1, tchr2, tstart2, tend2, inseqSize, tdir1, tdir2 = tkey.split('\t')

        if F[0] != tchr1 or int(F[1]) > int(tend1) + 25:

            print key2info[tkey]
            delList.append(tkey)
    
        else:

            if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2 and abs(int(F[2]) - len(tend1)) <= 25 and abs(int(F[5]) - len(tend2)) <= 25:

                infos = key2info[tkey].split('\t')
                if len(F[6].split(';')) < len(infos[6].split(';')) and len(F[6].split(';')) <= 2:
                    skipFlag = 1
                elif len(F[6].split(';')) > len(infos[6].split(';')) and len(infos[6].split(';')) <= 2:
                    delList.append(tkey)


    for tkey in delList:
        del key2info[tkey]

    if skipFlag == 0:
        key2info['\t'.join(F[0:6] + F[7:10])] = '\t'.join(F)

hIN.close()


for tkey in key2info:

    print key2info[tkey]

