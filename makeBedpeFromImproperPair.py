#! /usr/local/bin/python

"""
    script or creating bedpe file for representing possible breakpoint
"""

import sys

inputFile = sys.argv[1]

junction_dist = 500;
tempID, tempPairNum, tempChr, tempStart, tempEnd, tempDir, tempMapQ = "", "", "", "", "", "", ""

hIN = open(inputFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')

    pairNum = F[0][-1:]
    F[0] = F[0][0:-2]

    if F[0] == tempID and tempPairNum == "1" and pairNum == "2":
        
        chr1, dir1, start1, end1, mapQ1, align1 = tempChr, tempDir, 0, 0, tempMapQ, tempChr + ":" + str(tempStart) + "-" + str(tempEnd) 
        if dir1 == "+":
            start1 = tempEnd
            end1 = tempEnd + junction_dist
        else:
            start1 = tempStart - junction_dist
            end1 = tempStart

        chr2, dir2, start2, end2, mapQ2, align2 = F[1], F[4], 0, 0, F[5], F[1] + ":" + F[2] + "-" + F[3]
        if dir2 == "+":
            start2 = int(F[3])
            end2 = int(F[3]) + junction_dist
        else:
            start2 = int(F[2]) - junction_dist
            end2 = int(F[2])
 
        if chr1 < chr2:
            print '\t'.join([chr1, str(start1), str(end1), chr2, str(start2), str(end2), tempID, mapQ1 + "," + mapQ2, dir1, dir2, align1 + "," + align2])
        elif chr1 > chr2:
            print '\t'.join([chr2, str(start2), str(end2), chr1, str(start1), str(end1), tempID, mapQ2 + "," + mapQ1, dir2, dir1, align2 + "," + align1])
        else:
            if start1 <= start2:
                print '\t'.join([chr1, str(start1), str(end1), chr2, str(start2), str(end2), tempID, mapQ1 + "," + mapQ2, dir1, dir2, align1 + "," + align2])
            else:
                print '\t'.join([chr2, str(start2), str(end2), chr1, str(start1), str(end1), tempID, mapQ2 + "," + mapQ1, dir2, dir1, align2 + "," + align1])


    tempID, tempPairNum, tempChr, tempStart, tempEnd, tempDir, tempMapQ = F[0], pairNum, F[1], int(F[2]), int(F[3]), F[4], F[5]

hIN.close()

