#! /usr/local/bin/python

import sys, numpy, math
from scipy import stats


def summarizeRefAlt(inputFile):
 
    hIN = open(inputFile, 'r')

    numOther = 0
    numAlt = 0
    numRef = 0

    tempID = ""
    tempAlt = []
    tempRef1 = []
    tempRef2 = []
    tempRef = []
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        if F[0].isdigit() == False: continue

        F[9] = F[9][0:-2]
        if tempID != F[9]:
            if tempID != "":

                tempAltNM = 100
                if len(tempAlt) == 2:
                    if (tempAlt[0][1] <= tempAlt[1][1] and tempAlt[0][2] == "+" and tempAlt[1][2] == "-") or (tempAlt[0][1] >= tempAlt[1][1] and tempAlt[0][2] == "-" and tempAlt[1][2] == "+"): 
                        tempAltNM = tempAlt[0][0] + tempAlt[1][0]

                tempRefNM = 100
                if len(tempRef1) == 2:
                    if (tempRef1[0][1] <= tempRef1[1][1] and tempRef1[0][2] == "+" and tempRef1[1][2] == "-") or (tempRef1[0][1] >= tempRef1[1][1] and tempRef1[0][2] == "-" and tempRef1[1][2] == "+"): 
                        tempRefNM = min(tempRefNM, tempRef1[0][0] + tempRef1[1][0])

                if len(tempRef2) == 2:
                    if (tempRef2[0][1] <= tempRef2[1][1] and tempRef2[0][2] == "+" and tempRef2[1][2] == "-") or (tempRef2[0][1] >= tempRef2[1][1] and tempRef2[0][2] == "-" and tempRef2[1][2] == "+"): 
                        tempRefNM = min(tempRefNM, tempRef2[0][0] + tempRef2[1][0])

                if len(tempRef) == 2:
                    if (tempRef[0][1] <= tempRef[1][1] and tempRef[0][2] == "+" and tempRef[1][2] == "-") or (tempRef[0][1] >= tempRef[1][1] and tempRef[0][2] == "-" and tempRef[1][2] == "+"): 
                        tempRefNM = min(tempRefNM, tempRef[0][0] + tempRef[1][0])

                if tempAltNM >= 30 and tempRefNM >= 30:
                    numOther = numOther + 1
                elif tempAltNM < tempRefNM:
                    numAlt = numAlt + 1
                else:
                    numRef = numRef + 1


            tempID = F[9] 
            tempAlt = []
            tempRef1 = []
            tempRef2 = []
            tempRef = []

        tNM = int(F[10]) - int(F[0])
        tpos = int(F[15])
        tdir = F[8]

        if F[13][-3:] == "alt":
            tempAlt.append((tNM, tpos, tdir))
        elif F[13][-4:] == "ref1":
            tempRef1.append((tNM, tpos, tdir))
        elif F[13][-4:] == "ref2":
            tempRef2.append((tNM, tpos, tdir))
        elif F[13][-3:] == "ref":
            tempRef.append((tNM, tpos, tdir))


    tempAltNM = 100
    if len(tempAlt) == 2:
        if (tempAlt[0][1] <= tempAlt[1][1] and tempAlt[0][2] == "+" and tempAlt[1][2] == "-") or (tempAlt[0][1] >= tempAlt[1][1] and tempAlt[0][2] == "-" and tempAlt[1][2] == "+"):
            tempAltNM = tempAlt[0][0] + tempAlt[1][0]

    tempRefNM = 100
    if len(tempRef1) == 2:
        if (tempRef1[0][1] <= tempRef1[1][1] and tempRef1[0][2] == "+" and tempRef1[1][2] == "-") or (tempRef1[0][1] >= tempRef1[1][1] and tempRef1[0][2] == "-" and tempRef1[1][2] == "+"):
            tempRefNM = min(tempRefNM, tempRef1[0][0] + tempRef1[1][0])

    if len(tempRef2) == 2:
        if (tempRef2[0][1] <= tempRef2[1][1] and tempRef2[0][2] == "+" and tempRef2[1][2] == "-") or (tempRef2[0][1] >= tempRef2[1][1] and tempRef2[0][2] == "-" and tempRef2[1][2] == "+"):
            tempRefNM = min(tempRefNM, tempRef2[0][0] + tempRef2[1][0])

    if len(tempRef) == 2:
        if (tempRef[0][1] <= tempRef[1][1] and tempRef[0][2] == "+" and tempRef[1][2] == "-") or (tempRef[0][1] >= tempRef[1][1] and tempRef[0][2] == "-" and tempRef[1][2] == "+"):
            tempRefNM = min(tempRefNM, tempRef[0][0] + tempRef[1][0])


    if tempAltNM >= 30 and tempRefNM >= 30:
        numOther = numOther + 1
    elif tempAltNM < tempRefNM:
        numAlt = numAlt + 1
    else:
        numRef = numRef + 1

    return([numRef, numAlt])


tumorPsl = sys.argv[1]
normalPsl = sys.argv[2]

tumorRef, tumorAlt = summarizeRefAlt(tumorPsl)
normalRef, normalAlt = summarizeRefAlt(normalPsl)

# fisher test
oddsratio, pvalue = stats.fisher_exact([[tumorRef, tumorAlt], [normalRef, normalAlt]], 'less')

if pvalue < 1e-100:
    pvalue = 1e-100

lpvalue = (- math.log(pvalue, 10) if pvalue < 1 else 0)

print '\t'.join([str(tumorRef), str(tumorAlt), str(normalRef), str(normalAlt), str(lpvalue)])

 
