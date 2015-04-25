#! /usr/local/bin/python

import sys, numpy, math
from scipy import stats

"""
    note:
    Now the definition of "reference-prefer" read pair is those align at least 5 base than to alternative sequence. 
    but current definition generates problem when detecting mid-range tandem duplication on repeat sequences (5,115279667,-,5,115280072 (ATL-15T))
    because in this case, read pairs that should prefer the reference sequence do not prefer it significantly to alternative base...
    One remedy for this may be to 
    1. in advance, we should remove the read pairs whose source (reference or altenative) are unclear
    2. we define "class 2 reference read pair", and use them with "class 1 reference read pair".


"""


def checkScore(align):

    tempScore = 100
    if len(align) >= 2:
        for i1 in range(0, len(align) - 1):
            for i2 in range(i1 + 1, len(align)):

                if align[i1][1] <= align[i2][1] and align[i1][2] == "+" and align[i2][2] == "-":
                    tempScore = min(tempScore, align[i1][0] + align[i2][0])

                if align[i2][1] <= align[i1][1] and align[i2][2] == "+" and align[i1][2] == "-":
                    tempScore = min(tempScore, align[i1][0] + align[i2][0])

    return(tempScore)


def summarizeRefAlt(inputFile):
 
    hIN = open(inputFile, 'r')

    numOther = 0
    numAlt = 0
    numRef = 0

    # ref_ID = []
    # alt_ID = []
    # other_ID = []

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

                if tempID == "HWI-ST1021:128:C1K77ACXX:6:2303:6851:109312":
                    pass

                tempAltNM = checkScore(tempAlt)
                tempRefNM = min(checkScore(tempRef1), checkScore(tempRef2), checkScore(tempRef))

                if tempAltNM >= 30 and tempRefNM >= 30:
                    numOther = numOther + 1
                    # other_ID.append(tempID)
                elif tempAltNM < tempRefNM - 5:
                    numAlt = numAlt + 1
                    # alt_ID.append(tempID)
                elif tempRefNM < tempAltNM - 5:
                    numRef = numRef + 1
                    # ref_ID.append(tempID)


            tempID = F[9] 
            tempAlt = []
            tempRef1 = []
            tempRef2 = []
            tempRef = []

        tNM = int(F[10]) - int(F[0]) + int(F[5]) + int(F[7])
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


    tempAltNM = checkScore(tempAlt)
    tempRefNM = min(checkScore(tempRef1), checkScore(tempRef2), checkScore(tempRef))

    if tempAltNM >= 30 and tempRefNM >= 30:
        numOther = numOther + 1
        # other_ID.append(tempID)
    elif tempAltNM < tempRefNM - 5:
        numAlt = numAlt + 1
        # alt_ID.append(tempID)
    elif tempRefNM < tempAltNM - 5:
        numRef = numRef + 1
        # ref_ID.append(tempID)

    """
    print "ref"
    print '\n'.join(ref_ID)
    print "alt"
    print '\n'.join(alt_ID)
    print "other"
    print '\n'.join(other_ID)
    """

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

 
