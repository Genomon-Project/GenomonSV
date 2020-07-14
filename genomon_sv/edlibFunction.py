#!/usr/bin/env python

from __future__ import print_function
import sys, pysam
from . import utils
import edlib

def getRefAltForSV(inputFile):

    code = ""
    fa_alt = ""
    fa_ref1 = ""
    fa_ref2 = "" 
    fa_ref = ""
    with open(inputFile, "r") as hin:
        for line in hin:
            line = line.rstrip('\n')
            # Handling the eader 
            if line.startswith(">"):
                if line.endswith("alt"):
                    code = "alt"
                elif line.endswith("ref1"):
                    code = "ref1"
                elif line.endswith("ref2"):
                    code = "ref2"
                elif line.endswith("ref"):
                    code = "ref"
            # Handling the body 
            else:
                if code == "alt":
                    fa_alt = line
                elif code == "ref1":
                    fa_ref1 = line
                elif code == "ref2":
                    fa_ref2 = line
                elif code == "ref":
                    fa_ref = line
        
        return([fa_alt, fa_ref1, fa_ref2, fa_ref])
           

def summarizeRefAlt(inputFile, STDFlag, fa_alt, fa_ref1, fa_ref2, fa_ref, outputFile):

    read_pair = 0 
    name = "" 
    d_edlib_read1 = {}
    d_edlib_read2 = {}
    hout = open(outputFile, "w")
    with open(inputFile, "r") as hinT:
        for line in hinT:
            line = line.rstrip('\n')
            # Handling the eader 
            if line.startswith(">"):
                if line.endswith("/1"):
                    read_pair = 1 
                elif line.endswith("/2"):
                    read_pair = 2 
                else:
                    utils.warningMessage("ID in FASTA files do not have a /1 or /2 suffix. " + line)
                line = line.replace(">","",1).replace("/1","").replace("/2","")
                name = line  
            # Handling the body 
            else:
                rc_line = utils.reverseComplement(line)
                # edlib_ret
                # [0]=alt, [1]=alt_rev,
                # [2]=ref1,[3]=ref1_rev,
                # [4]=ref2,[5]=ref2_rev,
                # [6]=ref, [7]=ref_rev
                edlib_ret = [100,100,100,100,100,100,100,100]

                # alt
                ret = edlib.align(line, fa_alt, mode="HW", task="path")
                edlib_ret[0] = ret["editDistance"]
                ret = edlib.align(rc_line, fa_alt, mode="HW", task="path")
                edlib_ret[1] = ret["editDistance"]

                if fa_ref == "":
                    # ref1
                    ret = edlib.align(line, fa_ref1, mode="HW", task="path")
                    edlib_ret[2] = ret["editDistance"]
                    ret = edlib.align(rc_line, fa_ref1, mode="HW", task="path")
                    edlib_ret[3] = ret["editDistance"]
                    # ref2
                    ret = edlib.align(line, fa_ref2, mode="HW", task="path")
                    edlib_ret[4] = ret["editDistance"]
                    ret = edlib.align(rc_line, fa_ref2, mode="HW", task="path")
                    edlib_ret[5] = ret["editDistance"]
                else:
                    # ref
                    ret = edlib.align(line, fa_ref, mode="HW", task="path")
                    edlib_ret[6] = ret["editDistance"]
                    ret = edlib.align(rc_line, fa_ref, mode="HW", task="path")
                    edlib_ret[7] = ret["editDistance"]

                if read_pair == 1:
                    d_edlib_read1[name] = edlib_ret
                else:
                    d_edlib_read2[name] = edlib_ret

    numRef = 0 
    numAlt = 0
    numOther = 0
    for name in d_edlib_read1.keys():
        ret1 = d_edlib_read1[name]
        ret2 = d_edlib_read2[name]

        ed_alt_tmp_1  = ret1[0]+ret2[1]
        ed_alt_tmp_2  = ret1[1]+ret2[0]
        ed_alt = min(ed_alt_tmp_1, ed_alt_tmp_2)

        ed_ref = 0
        if fa_ref == "":
            ed_ref1_tmp_1 = ret1[2]+ret2[3]
            ed_ref1_tmp_2 = ret1[3]+ret2[2]
            ed_ref2_tmp_1 = ret1[4]+ret2[5]
            ed_ref2_tmp_2 = ret1[5]+ret2[4]
            ed_ref = min(ed_ref1_tmp_1, ed_ref1_tmp_2, ed_ref2_tmp_1, ed_ref2_tmp_2)
            # print(name ,file=hout)
            # print("\t".join([str(ed_alt_tmp_1), str(ed_alt_tmp_2), str(ed_ref1_tmp_1), str(ed_ref1_tmp_2), str(ed_ref2_tmp_1), str(ed_ref2_tmp_2)]), file=hout)
        else:
            ed_ref_tmp_1  = ret1[6]+ret2[7]
            ed_ref_tmp_2  = ret1[7]+ret2[6]
            ed_ref = min(ed_ref_tmp_1, ed_ref_tmp_2)
            # print(name ,file=hout)
            # print("\t".join([str(ed_alt_tmp_1), str(ed_alt_tmp_2), str(ed_ref_tmp_1), str(ed_ref_tmp_2)]), file=hout)

        # ed_alt  = min(ret1[0],ret1[1]) + min(ret2[0],ret2[1])
        # ed_ref1_tmp = min(ret1[2],ret1[3]) + min(ret2[2],ret2[3])
        # ed_ref2_tmp = min(ret1[4],ret1[5]) + min(ret2[4],ret2[5])
        # ed_ref_tmp  = min(ret1[6],ret1[7]) + min(ret2[6],ret2[7])
        # ed_ref = min(ed_ref1_tmp, ed_ref2_tmp, ed_ref_tmp)

        if ed_alt >= 10 and ed_ref >= 10:
            numOther += 1
        elif ed_alt < ed_ref - 5:
            numAlt += 1
        elif STDFlag == 0 and ed_ref < ed_alt - 5 or STDFlag == 1 and ed_ref <= ed_alt:
            numRef += 1 

    hout.close()
    return([numRef, numAlt])


