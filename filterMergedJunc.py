#! /usr/local/bin/python

# script for filtering the inferred break points for structural variation
#
# as conditions for filtering
# the number of support reads
# mapping quality (e.g., in either type 1 or 2 junction reads, median of mapping quality is at least 40 ??)
# the size of alignment covered region (e.g., around the both break points, at least 80 bp sequences have to be covered)

import sys
import numpy
import re

inputFile = sys.argv[1]
minSupReadNum = int(sys.argv[2])
minMapQ = int(sys.argv[3])
minCoverBase = int(sys.argv[4])

##########
# functions and classes for comparing regions
regionRe = re.compile('(\w+):(\d+)\-(\d+)')

class Regions(object):

    regionVec = []
    def __init__(self):
        self.regionVec = []

    def addMerge(self, newRegion):
        mergedFlag = 0    
        for i in range(0, len(self.regionVec)):
            tempReg = self.regionVec[i]
            newReg = regionMerge(self.regionVec[i], newRegion)
            if newReg != '':
                mergedFlag = 1
                self.regionVec[i] = newReg

        if mergedFlag == 0:
            self.regionVec.append(newRegion)    

    def checkMerge(self, newRegion):
        mergedFlag = 0
        for i in range(0, len(self.regionVec)):
            tempReg = self.regionVec[i]
            newReg = regionMerge(self.regionVec[i], newRegion)
            if newReg != '':
                mergedFlag = 1
                self.regionVec[i] = newReg

    def reduceMerge(self):
        mergedFlag = 1
        while mergedFlag == 1:
            mergedFlag = 0
            mergePair = []
            for i in range(0, len(self.regionVec) - 1):
                for j in range(i + 1, len(self.regionVec)):
                    newReg = regionMerge(self.regionVec[i], self.regionVec[j])
                    if newReg != '':
                        mergedFlag = 1 
                        mergePair = [i, j]
                        
            if mergedFlag == 1:
                self.regionVec[mergePair[0]] = regionMerge(self.regionVec[mergePair[0]], self.regionVec[mergePair[1]])
                del self.regionVec[mergePair[1]]



    def regionSize(self):
        size = 0
        for reg in self.regionVec:
            reMatch = regionRe.match(reg)
            chr = reMatch.group(1)
            start = int(reMatch.group(2))
            end = int(reMatch.group(3))
            size += end - start + 1
        return(size)


def regionMerge(region1, region2):

    reMatch1 = regionRe.match(region1)
    reMatch2 = regionRe.match(region2)
    chr1 = reMatch1.group(1)
    start1 = int(reMatch1.group(2))
    end1 = int(reMatch1.group(3))
    chr2 = reMatch2.group(1)
    start2 = int(reMatch2.group(2))
    end2 = int(reMatch2.group(3))
    
    # if there is no overlap, then return region for the 
    if chr1 != chr2 or end1 < start2 or end2 < start1:
        return '' 
        
    return chr1 + ":" + str(min(start1, start2)) + "-" + str(max(end1, end2))
####################
##########

inFile = open(inputFile, 'r')


for line in inFile:
    F = line.rstrip('\n').split('\t')

    ##########
    # for now only consider long ranged SV??
    if F[0] == F[3] and abs(int(F[2]) - int(F[5])) < 1000:
        continue
    ##########


    MQs = F[7].split(';')
    coveredRegions = F[10].split(';')
    pairCoveredRegions = F[13].split(';')
    pairMQs = F[12].split(';')
    juncReadTypes = F[14].split(';')
    improperCoveredRegion = ([] if F[17] == "---" else F[17].split(';'))

    # enumerate support read number
    juncSupport = len(MQs)
    improperSupport = (0 if F[16] == "---" else len(F[16].split(';')))
    # skip if the number of suppor read is below the minSupReadNum 
    if juncSupport + improperSupport < minSupReadNum:
        continue

    ##########
    # check for the mapping quality
    MQs1 = [];
    MQs2 = [];
    for i in range(0, len(MQs)):
        if juncReadTypes[i] == "1":
            MQs1.append(int(MQs[i]))
        else:
            MQs2.append(int(MQs[i]))

    
    mapFlag = 0
    if len(MQs1) > 0 and numpy.median(MQs1) >= minMapQ:
        mapFlag = 1

    if len(MQs2) > 0 and numpy.median(MQs2) >= minMapQ:
        mapFlag = 1     

    if len(pairMQs) > 0 and numpy.median(pairMQs) >= minMapQ:
        mapFlag = 1

    if mapFlag == 0:
        continue
    ##########


    ########## 
    # check for the covered region
    region1 = Regions()
    region2 = Regions()

    for i in range(0, len(coveredRegions)):
        regPair = coveredRegions[i].split(',')
        if juncReadTypes[i] == "1":
            region1.addMerge(regPair[0])
            region2.addMerge(regPair[1])
        else:
            region1.addMerge(regPair[1])
            region2.addMerge(regPair[0])

    targetRegion1 = ""
    if F[8] == "+":
        targetRegion1 = F[0] + ":" + str(max(0, int(F[2]) - 1000)) + "-" + F[2]
    else:
        targetRegion1 = F[0] + ":" + F[2] + "-" + str(int(F[2]) + 1000)

    targetRegion2 = ""
    if F[9] == "+":
        targetRegion2 = F[3] + ":" + str(max(0, int(F[5]) - 1000)) + "-" + F[5]
    else:
        targetRegion2 = F[3] + ":" + F[5] + "-" + str(int(F[5]) + 1000)

    for i in range(0, len(pairCoveredRegions)):
        if regionMerge(targetRegion1, pairCoveredRegions[i]):
            region1.addMerge(pairCoveredRegions[i])        
        elif regionMerge(targetRegion2, pairCoveredRegions[i]):
            region2.addMerge(pairCoveredRegions[i])

    for i in range(0, len(improperCoveredRegion)):
        regPair = improperCoveredRegion[i].split(',')
        if regionMerge(targetRegion1, regPair[0]) and regionMerge(targetRegion2, regPair[1]):
            region1.addMerge(regPair[0])
            region2.addMerge(regPair[1])
        elif regionMerge(targetRegion1, regPair[1]) and regionMerge(targetRegion2, regPair[0]):
            region1.addMerge(regPair[1])
            region2.addMerge(regPair[0])

    # print >> sys.stderr, F[6] 
    region1.reduceMerge()
    region2.reduceMerge()

    if region1.regionSize() < minCoverBase or region2.regionSize() < minCoverBase:
        continue

    # every condition is satisfied
    print '\t'.join(F)


inFile.close()

