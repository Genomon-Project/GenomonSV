#! /usr/local/bin/python

"""
    script for obtaining pair read information (mainly end position, because it cannot recovered from bam files)
"""

import sys, re, pysam, tabix

inputBAM = sys.argv[1]
inputBedpe = sys.argv[2]
region = sys.argv[3]

bamfile = pysam.Samfile(inputBAM, "rb")
# tabixfile = pysam.Tabixfile(inputBedpe)
tabixfile = tabix.open(inputBedpe)
 
regionMatch = re.search(r'(\w+):(\d+)\-(\d+)', region)
chr_region = regionMatch.group(1)
start_region = regionMatch.group(2)
end_region = regionMatch.group(3)

ID2info = {}
# for row in tabixfile.fetch(chr_region, int(start_region), int(end_region)):
#     F = row.rstrip('\n').split('\t')
#     ID2info[F[3]] = '\t'.join(F)

tabixErrorFlag = 0
try:
    records = tabixfile.query(chr_region, int(start_region), int(end_region))
except Exception as inst:
    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
    tabixErrorFlag = 1

normalJunctionIDs = [];
if tabixErrorFlag == 0:
    for record in records:
        ID2info[record[3]] = '\t'.join(record)


for read in bamfile.fetch(chr_region, int(start_region), int(end_region)):

    flags = format(int(read.flag), '#014b')[:1:-1]

    # skip supplementary alignment
    if flags[8] == "1" or flags[11] == "1": continue

    # skip one of the pair is unmapped
    if flags[2] == "1" or flags[3] == "1": continue
 
    seqID = (read.qname + "/1" if  flags[6] == "1" else read.qname + "/2")

    # if re.match("HWI-ST1021:119:C14DVACXX:1:2301:9987:12915", read.qname):
        # print read.qname
 
    if seqID in ID2info:
        print ID2info[seqID] + "\t" + bamfile.getrname(read.tid) + ":" + str(read.pos + 1) + "-" + str(read.aend) + "\t" + str(read.mapq)



