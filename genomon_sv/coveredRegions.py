#!/usr/bin/env python

import re

regionRe = re.compile('([^ \t\n\r\f\v,]+):(\d+)\-(\d+)')


class coveredRegions(object):

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

