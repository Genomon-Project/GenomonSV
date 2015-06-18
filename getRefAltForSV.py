#! /usr/local/bin/python

"""
    for short SV (mid-range (<= 1000bp) deletion, tandem duplication), we get the two sequence
    for large SV (>1000bp), we get three sequence (one joint sequence by two break points, and two reference sequences around the break points)

    the only concern is short inversion... (are there some somatic short inversion?)
    however, this will be filtered beforehand by the "cover filter", and maybe we have to give up detecting these class of SVs.

"""

import pysam

def getRefAltForSV(reference, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, juncSeq):

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    if juncSeq == "---": juncSeq = ""

    # for mid-range deletion or tandemu duplication
    if juncChr1 == juncChr2 and abs(juncPos1 - juncPos2) <= 1000 and juncDir1 != juncDir2:

        seq = ""
        for item in pysam.faidx(reference, juncChr1 + ":" + str(juncPos1 - 1000) + "-" + str(juncPos2 + 1000)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()

        print '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref"
        print seq

        # for mid-range deletion
        if juncDir1 == "+" and juncDir2 == "-":

            seq = ""
            for item in pysam.faidx(reference, juncChr1 + ":" + str(juncPos1 - 1000) + "-" + str(juncPos1)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()

            seq = seq + juncSeq

            for item in pysam.faidx(reference, juncChr2 + ":" + str(juncPos2) + "-" + str(juncPos2 + 1000)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()

            print '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt"
            print seq

        # for mid-range tandem duplication
        else:
            seq = "" 
            for item in pysam.faidx(reference, juncChr2 + ":" + str(juncPos2 - 1000) + "-" + str(juncPos2)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()
            
            seq = seq + juncSeq

            for item in pysam.faidx(reference, juncChr1 + ":" + str(juncPos1) + "-" + str(juncPos1 + 1000)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()
            
            print '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt"
            print seq



    else:

        seq = ""
        for item in pysam.faidx(reference, juncChr1 + ":" + str(juncPos1 - 1000) + "-" + str(juncPos1 + 1000)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()

        print '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref1"
        print seq

        seq = ""
        for item in pysam.faidx(reference, juncChr2 + ":" + str(juncPos2 - 1000) + "-" + str(juncPos2 + 1000)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()
            
        print '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref2"
        print seq


       
        seq = ""
        if juncDir1 == "+":
            tseq = ""
            for item in pysam.faidx(reference, juncChr1 + ":" + str(juncPos1 - 1000) + "-" + str(juncPos1)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
        else:
            tseq = ""
            for item in pysam.faidx(reference, juncChr1 + ":" + str(juncPos1) + "-" + str(juncPos1 + 1000)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
            tseq = "".join(complement.get(base) for base in reversed(tseq))

        seq = tseq + juncSeq

        if juncDir2 == "-":
            tseq = "" 
            for item in pysam.faidx(reference, juncChr2 + ":" + str(juncPos2) + "-" + str(juncPos2 + 1000)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
        else:
            tseq = ""
            for item in pysam.faidx(reference, juncChr2 + ":" + str(juncPos2 - 1000) + "-" + str(juncPos2)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
            tseq = "".join(complement.get(base) for base in reversed(tseq))

        seq = seq + tseq

        print '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt"
        print seq



if __name__ == "__main__":
    import sys, pysam
    getRefAltForSV(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5], int(sys.argv[6]), sys.argv[7], sys.argv[8])


