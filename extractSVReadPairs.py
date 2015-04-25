#! /usr/local/bin/python


"""
    read pairs containing break points are extracted. (yshira 2015/04/23)
    The exact condition is as follows:

    1. one of the read in the pair has the break point of the SV candidate
    2. the start positions of the read pairs are within 800bp of the break point of the SV candidate 

    Some minor concern for the above conditions are:
    1. Depending on the choice of the "start position" or "end position", the distance between the read and break point differs. This can generate slight bias...
    (but I believe we can recover this by setting sufficient margin (800bp), and summarize the alignment result carefully.)
    2. Maybe, for some of the read pair, the result of alignment is obvious. But should we re-align them?

"""

import pysam


maxDepth = 1000
# get the read pairs covering the junction break point

def extractSVReadPairs(bamFilePath, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, searchLength, margin):

    bamfile = pysam.Samfile(bamFilePath, 'rb')

    # if the #sequence read is over the `maxDepth`, then that key is ignored
    depthFlag = 0
    if bamfile.count(juncChr1, int(juncPos1) - 1, int(juncPos1) + 1) >= maxDepth: depthFlag = 1
    if bamfile.count(juncChr2, int(juncPos2) - 1, int(juncPos2) + 1) >= maxDepth: depthFlag = 1
    if depthFlag == 1: 
        sys.exit(27)

    readIDs = []    
    for read in bamfile.fetch(juncChr1, int(juncPos1) - searchLength, int(juncPos1) + searchLength):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue 

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")
        chr_pair = bamfile.getrname(read.rnext)
        pos_pair = int(read.pnext + 1)
        dir_pair = ("-" if flags[5] == "1" else "+")

        # the read (with margin) contains break point
        if pos_current - margin <= juncPos1 <= (read.aend - 1) + margin:
            readIDs.append(read.qname)
    
        # the read pair covers break point
        if chr_pair == juncChr1 and pos_current <= juncPos1 <= pos_pair and dir_current == "+" and dir_pair == "-":
            readIDs.append(read.qname)

        # the read pair covers break point
        if chr_pair == juncChr2:
            juncFlag = 0
            if juncDir1 == "+" and juncDir2 == "+" and pos_current <= juncPos1 and pos_pair <= juncPos2: juncFlag = 1
            if juncDir1 == "+" and juncDir2 == "-" and pos_current <= juncPos1 and pos_pair >= juncPos2: juncFlag = 1
            if juncDir1 == "-" and juncDir2 == "+" and pos_current >= juncPos1 and pos_pair <= juncPos2: juncFlag = 1
            if juncDir1 == "-" and juncDir2 == "-" and pos_current >= juncPos1 and pos_pair >= juncPos2: juncFlag = 1

            if juncFlag == 1:  
               readIDs.append(read.qname)


    for read in bamfile.fetch(juncChr2, int(juncPos2) - searchLength, int(juncPos2) + searchLength):
        
        if read.qname == "ST-E00104:162:H03UUALXX:5:1222:21168:16006":
            pass
 
        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue
        
        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue
        
        # skip duplicated reads
        if flags[10] == "1": continue
        
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")
        chr_pair = bamfile.getrname(read.rnext)
        pos_pair = int(read.pnext + 1)
        dir_pair = ("-" if flags[5] == "1" else "+")
        
        # the read (with margin) contains break point
        if pos_current - margin <= juncPos2 <= (read.aend - 1) + margin:
            readIDs.append(read.qname)
                
        # the read pair covers break point
        if chr_pair == juncChr2 and pos_current <= juncPos2 <= pos_pair and dir_current == "+" and dir_pair == "-":
            readIDs.append(read.qname)
                
        # the read pair covers break point
        if chr_pair == juncChr1:
            juncFlag = 0
            if juncDir2 == "+" and juncDir1 == "+" and pos_current <= juncPos2 and pos_pair <= juncPos1: juncFlag = 1
            if juncDir2 == "+" and juncDir1 == "-" and pos_current <= juncPos2 and pos_pair >= juncPos1: juncFlag = 1
            if juncDir2 == "-" and juncDir1 == "+" and pos_current >= juncPos2 and pos_pair <= juncPos1: juncFlag = 1
            if juncDir2 == "-" and juncDir1 == "-" and pos_current >= juncPos2 and pos_pair >= juncPos1: juncFlag = 1
             
            if juncFlag == 1:
               readIDs.append(read.qname)


    readID2seq1 = {}
    readID2seq2 = {}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for read in bamfile.fetch(juncChr1, int(juncPos1) - searchLength, int(juncPos1) + searchLength):

        if read.qname in readIDs:
        
            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip unmapped read 
            if flags[2] == "1" or flags[3] == "1": continue

            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue

            # skip duplicated reads
            if flags[10] == "1": continue

            tempSeq = ""
            if flags[4] == "1":
                tempSeq = "".join(complement.get(base) for base in reversed(str(read.seq)))
            else:
                tempSeq = read.seq
 
            # the first read
            if flags[6] == "1":
                readID2seq1[read.qname] = tempSeq
            else:
                readID2seq2[read.qname] = tempSeq

                 
    for read in bamfile.fetch(juncChr2, int(juncPos2) - searchLength, int(juncPos2) + searchLength):

        if read.qname in readIDs:

            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip unmapped read 
            if flags[2] == "1" or flags[3] == "1": continue

            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue
            
            # skip duplicated reads
            if flags[10] == "1": continue

            tempSeq = ""
            if flags[4] == "1":
                tempSeq = "".join(complement.get(base) for base in reversed(str(read.seq)))
            else:
                tempSeq = read.seq

            # the first read
            if flags[6] == "1":
                readID2seq1[read.qname] = tempSeq
            else:
                readID2seq2[read.qname] = tempSeq


    for readID in readID2seq1:
        if readID in readID2seq2:
            print '>' + readID + '/1'
            print readID2seq1[readID]
            print '>' + readID + '/2'
            print readID2seq2[readID]



if __name__ == "__main__":
    import sys, pysam
    extractSVReadPairs(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5], int(sys.argv[6]), sys.argv[7], int(sys.argv[8]), int(sys.argv[9]))


