#!/usr/bin/env python

"""
    functions for parsing breakpoint containing read pairs and improperly aligned read pairs
"""

import sys, pysam, re, subprocess, collections
import utils

def parseJunctionFromBam(inputBAM, outputFilePath, min_mapping_qual, abnormal_insert_size, min_major_clip_size, max_minor_clip_size):

    """
    This function utilizes the SA tags (SA:Z:rname, pos, strand, CIGAR, mapQ, number of mismatch).
    The strand of the supplementary alignment information in the SA tag is determined by the orignal sequence (before taking complement).
    Therefore, please not that when the primary alignment is in the reverse direction, the sequence shown in the bam file does not match
    to the SA tags..
    """

    bamfile = pysam.Samfile(inputBAM, "rb")
    hOUT = open(outputFilePath, "w")
 
    SAre = re.compile('([^ \t\n\r\f\v,]+),(\d+),([\-\+]),(\w+),(\d+),(\d+);')
    cigarMDRe = re.compile('(\d+)([MD])')
    cigarHIMSRe = re.compile('(\d+)([HIMS])')
    cigarHSRe_right = re.compile('(\d+)([HS])$')
    cigarHSRe_left = re.compile('^(\d+)([HS])')


    # maybe add the regional extraction of bam files
    for read in bamfile.fetch():

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip if either of the read pair is unmapped
        if flags[2] == "1" or flags[3] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        # no clipping
        if len(read.cigar) == 1: continue

        # skip if below the minimum mapping quality
        if (read.mapq < min_mapping_qual): continue

        # skip if the read aligned to hs37d5"
        # (in the future, this step will be replaced to some more sophisticated way;
        # (e.g., the user can input the target chromosomes and ignore if the read is aligned to non-target chromosomes, and so on..
        if bamfile.getrname(read.tid) == "hs37d5" or bamfile.getrname(read.rnext) == "hs37d5": continue


        # get the clipping size in the both side
        left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
        right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)

        if left_clipping < min_major_clip_size and right_clipping < min_major_clip_size: continue

        # for comparing with the previous script (this will removed soon)
        if left_clipping > max_minor_clip_size and right_clipping > max_minor_clip_size: continue

        # skip if there is no SA tags
        SA_str = None 
        for item in read.tags:
            if item[0] == "SA":
                SA_str = SAre.match(item[1])
                
        if SA_str is None: continue

        # get the alignment basic information
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")
        chr_pair = bamfile.getrname(read.rnext)
        pos_pair = int(read.pnext + 1)
        dir_pair = ("-" if flags[5] == "1" else "+")

        chr_SA, pos_SA, dir_SA, cigar_SA = SA_str.group(1), int(SA_str.group(2)), SA_str.group(3), SA_str.group(4)


        # get the soft clipping information on the supplementary alignment
        right_clipping_SA = 0
        tempMatch = cigarHSRe_right.search(cigar_SA)
        if tempMatch is not None: right_clipping_SA = int(tempMatch.group(1))

        left_clipping_SA = 0
        tempMatch = cigarHSRe_left.search(cigar_SA)
        if tempMatch is not None: left_clipping_SA = int(tempMatch.group(1))

        # skip if the both sides of the supplementary alignemt is clipped than the specified threshould
        if left_clipping_SA > max_minor_clip_size and right_clipping_SA > max_minor_clip_size: continue


        # when the right side is clipped...
        if right_clipping >= min_major_clip_size:

            clipLen_current = right_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen

            juncChr_current = chr_current
            juncPos_current = pos_current + alignmentSize_current - 1
            juncDir_current = "+"
            coverRegion_current = chr_current + ":" + str(pos_current) + "-" + str(juncPos_current)
            juncChr_SA = chr_SA

            # get the expected clipped size in the supplementary read and the side of the clipping
            # for getting the expected side of the clipping, we have to be carefully consult with the definition of SA tag.
            expected_clipLen_SA = readLength_current - clipLen_current
            expected_clipDir_SA = ("-" if dir_current == dir_SA else "+")

            alignmentSize_SA = 0
            for item in cigarMDRe.finditer(cigar_SA):
                alignmentSize_SA += int(item.group(1))        


            validFlag = 0
            juncDir_SA = ""
            juncPos_SA = ""
            # the pair read is aligned at the same chromosome with the primary read
            if dir_current == "-" and dir_pair == "+" and chr_current == chr_pair and 0 <= pos_current - pos_pair < abnormal_insert_size:

                if dir_SA == "+" and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                    clipLen_SA = right_clipping_SA
                    juncDir_SA = "+"
                    juncPos_SA = pos_SA + alignmentSize_SA - 1
                    if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA) 
                    coverRegion_SA = chr_SA + ":" + str(pos_SA) + "-" + str(juncPos_SA)
                    juncType = 1
                    validFlag = 1

                if dir_SA == "-" and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                    clipLen_SA = left_clipping_SA 
                    juncDir_SA = "-"
                    juncPos_SA = pos_SA
                    if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                    coverRegion_SA = chr_SA + ":" + str(juncPos_SA) + "-" + str(juncPos_SA + alignmentSize_SA - 1)
                    juncType = 1
                    validFlag = 1

            # when the supplementary read is aligned on the same chromosome with the paired read
            if dir_current == "+" and chr_SA == chr_pair:

                if dir_SA == "-" and dir_pair == "+" and 0 <= pos_SA - pos_pair < abnormal_insert_size and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                    clipLen_SA = right_clipping_SA
                    juncDir_SA = "+"
                    juncPos_SA = pos_SA + alignmentSize_SA - 1
                    if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                    coverRegion_SA = chr_SA + ":" + str(pos_SA) + "-" + str(juncPos_SA) 
                    juncType = 2 
                    validFlag = 1

                if dir_SA == "+" and dir_pair == "-" and 0 <= pos_pair - pos_SA < abnormal_insert_size and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                    clipLen_SA = left_clipping_SA
                    juncDir_SA = "-"
                    juncPos_SA = pos_SA
                    if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                    coverRegion_SA = chr_SA + ":" + str(juncPos_SA) + "-" + str(juncPos_SA + alignmentSize_SA - 1)
                    juncType = 2
                    validFlag = 1


            if validFlag == 1:

                if juncPos_SA < 0: continue

                juncSurplus = "---"
                if clipLen_SA > expected_clipLen_SA and readLength_current == len(read.seq):
                    surPlus_start = readLength_current - clipLen_current
                    surPlus_end = surPlus_start + clipLen_SA - expected_clipLen_SA
                    juncSurplus = read.seq[surPlus_start:surPlus_end]

                # reorder by the chromosome position and print
                if juncChr_current < juncChr_SA or juncChr_current == juncChr_SA and juncPos_current <= juncPos_SA:
                    print >> hOUT, '\t'.join([juncChr_current, str(juncPos_current - 1), str(juncPos_current), juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), \
                                     read.qname + ("/1" if flags[6] == "1" else "/2"), juncSurplus, juncDir_current, juncDir_SA, \
                                     str(read.mapq), coverRegion_current + "," + coverRegion_SA, chr_pair + ":" + str(pos_pair), str(juncType), "1"])
                else: 
                    print >> hOUT, '\t'.join([juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), juncChr_current, str(juncPos_current - 1), str(juncPos_current), \
                                     read.qname + ("/1" if flags[6] == "1" else "/2"), juncSurplus, juncDir_SA, juncDir_current, \
                                     str(read.mapq), coverRegion_current + "," + coverRegion_SA, chr_pair + ":" + str(pos_pair), str(juncType), "2"])


        if left_clipping >= min_major_clip_size:

            clipLen_current = left_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen
     
            juncChr_current = chr_current
            juncPos_current = pos_current
            juncDir_current = "-"
            coverRegion_current = chr_current + ":" + str(pos_current) + "-" + str(pos_current + alignmentSize_current - 1)
            juncChr_SA = chr_SA

            expected_clipLen_SA = readLength_current - clipLen_current
            expected_clipDir_SA = ("+" if dir_current == dir_SA else "-")

            alignmentSize_SA = 0
            for item in cigarMDRe.finditer(cigar_SA):
                alignmentSize_SA += int(item.group(1))


            validFlag = 0
            juncDir_SA = ""
            juncPos_SA = ""
            # the pair read is aligned at the same chromosome with the primary read
            if dir_current == "+" and dir_pair == "-" and chr_current == chr_pair and 0 <= pos_pair - pos_current < abnormal_insert_size:

                if dir_SA == "+" and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                    clipLen_SA = right_clipping_SA
                    juncDir_SA = "+"
                    juncPos_SA = pos_SA + alignmentSize_SA - 1
                    if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                    coverRegion_SA = chr_SA + ":" + str(pos_SA) + "-" + str(juncPos_SA)
                    juncType = 1
                    validFlag = 1

                if dir_SA == "-" and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                    clipLen_SA = left_clipping_SA
                    juncDir_SA = "-"
                    juncPos_SA = pos_SA
                    if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                    coverRegion_SA = chr_SA + ":" + str(juncPos_SA) + "-" + str(juncPos_SA + alignmentSize_SA - 1)
                    juncType = 1
                    validFlag = 1


            if dir_current == "-" and chr_SA == chr_pair:

               if dir_SA == "-" and dir_pair == "+" and 0 <= pos_SA - pos_pair < abnormal_insert_size and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                    clipLen_SA = right_clipping_SA
                    juncDir_SA = "+"
                    juncPos_SA = pos_SA + alignmentSize_SA - 1
                    if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                    coverRegion_SA = chr_SA + ":" + str(pos_SA) + "-" + str(juncPos_SA)
                    juncType = 2
                    validFlag = 1

               if dir_SA == "+" and dir_pair == "-" and 0 <= pos_pair - pos_SA < abnormal_insert_size and expected_clipDir_SA == "-" and left_clipping_SA:
                    clipLen_SA = left_clipping_SA
                    juncDir_SA = "-"
                    juncPos_SA = pos_SA
                    if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                    coverRegion_SA = chr_SA + ":" + str(juncPos_SA) + "-" + str(juncPos_SA + alignmentSize_SA - 1)
                    juncType = 2
                    validFlag = 1


            if validFlag == 1:

                if juncPos_SA < 0: continue

                juncSurplus = "---"
                if clipLen_SA > expected_clipLen_SA and readLength_current == len(read.seq):
                    # surPlus_end = clipLen_current - 1 # this is not correct. but there seems to be a bug in the original perl script
                    surPlus_end = clipLen_current # this is right
                    surPlus_start = surPlus_end - (clipLen_SA - expected_clipLen_SA)
                    juncSurplus = read.seq[surPlus_start:surPlus_end]

                # reorder by the chromosome position and print
                if juncChr_current < juncChr_SA or juncChr_current == juncChr_SA and juncPos_current <= juncPos_SA: 
                    print >> hOUT, '\t'.join([juncChr_current, str(juncPos_current - 1), str(juncPos_current), juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), \
                                     read.qname + ("/1" if flags[6] == "1" else "/2"), juncSurplus, juncDir_current, juncDir_SA, \
                                     str(read.mapq), coverRegion_current + "," + coverRegion_SA, chr_pair + ":" + str(pos_pair), str(juncType), "1"])
                else:                
                    print >> hOUT, '\t'.join([juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), juncChr_current, str(juncPos_current - 1), str(juncPos_current), \
                                     read.qname + ("/1" if flags[6] == "1" else "/2"), juncSurplus, juncDir_SA, juncDir_current, \
                                     str(read.mapq), coverRegion_current + "," + coverRegion_SA, chr_pair + ":" + str(pos_pair), str(juncType), "2"])

    bamfile.close()
    hOUT.close()



def getPairStartPos(inputFilePath, outputFilePath):

    """
        script for obtaining the position information about the pair read from the junction file 

    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath + ".tmp", 'w')
    num = 1

    reChrPos = re.compile('^([^ \t\n\r\f\v,]+):(\d+)')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        ID = F[6]
        chr = ""
        pos = ""

        # obtain the information about the start site of the pair read 
        chrpos = reChrPos.search(F[12])
        if chrpos is not None:
            chr = chrpos.group(1)
            pos = chrpos.group(2)
        else:
            print '\t'.join(F)
            print "the 13th column did not match to (chr):(start) pattern"
            sys.exit()

        # change the pair read num
        if ID[-1:] == "1":
            ID = ID[:-1] + "2"
        else:
            ID = ID[:-1] + "1"

        print >> hOUT, chr + '\t' + str(int(pos) - 1) + '\t' + pos + '\t' + ID + '\t' + str(num)

        num = num + 1

    hIN.close()
    hOUT.close()


    hOUT = open(outputFilePath, 'w')
    subprocess.call(["sort", "-k1,1", "-k3,3n", outputFilePath + ".tmp"], stdout = hOUT)
    hOUT.close()


    ####################
    # delete intermediate file
    subprocess.call(["rm", outputFilePath + '.tmp'])



def getPairCoverRegionFromBam(inputBam, outputFilePath, inputTabixFile):

    """
        script for obtaining pair read information (mainly end position, because it cannot recovered from bam files)
    """
    ####################
    bamfile = pysam.Samfile(inputBam, "rb")
    tabixfile = pysam.TabixFile(inputTabixFile)
    hOUT = open(outputFilePath + ".tmp", "w")

    ID2info = {}
    tempChr = ""
    tempPos = 0
    checkPositionMargin = 10000000

    tabixErrorMsg = ""
    for read in bamfile.fetch():

        # when into new regions, fetch the keys from the tabix indexed file
        if bamfile.getrname(read.tid) != tempChr or int(read.pos + 1) > tempPos + checkPositionMargin:

            tempChr = bamfile.getrname(read.tid)
            tempPos = int(read.pos + 1) - 1

            ID2info = {}
            tabixErrorFlag = 0
            try:
                records = tabixfile.fetch(tempChr, tempPos, tempPos + checkPositionMargin)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorMsg = str(inst.args) 
                tabixErrorFlag = 1

            if tabixErrorFlag == 0:
                for record in records:
                    splt_record = record.split('\t')
                    ID2info[splt_record[3]] = record


        flags = format(int(read.flag), '#014b')[:1:-1]

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip one of the pair is unmapped
        if flags[2] == "1" or flags[3] == "1": continue
     
        seqID = (read.qname + "/1" if  flags[6] == "1" else read.qname + "/2")

     
        if seqID in ID2info:
            print >> hOUT, ID2info[seqID] + "\t" + bamfile.getrname(read.tid) + ":" + str(read.pos + 1) + "-" + str(read.aend) + "\t" + str(read.mapq)


    if tabixErrorMsg != "":
        utils.warningMessage("One or more error occured in tabix file fetch, e.g.: " + tabixErrorMsg)

    bamfile.close()
    tabixfile.close()
    ####################

    ####################
    hOUT = open(outputFilePath, 'w')
    subprocess.call(["sort", "-k5n", outputFilePath + ".tmp"], stdout = hOUT)
    hOUT.close()
    ####################


    ####################
    subprocess.call(["rm", outputFilePath + ".tmp"])



def addPairCoverRegionFromBam(inputFilePath, outputFilePath, pairCoverRegionInfoFile):

    """
        script for adding junction pair position and mapping quality information
    """

    hOriginalFile = open(inputFilePath, 'r')
    hPairInfoFile = open(pairCoverRegionInfoFile, 'r')
    hOUT = open(outputFilePath, 'w')

    line1 = hOriginalFile.readline().rstrip('\n')
    # if no line in the input file
    if line1  == '':
        hOriginalFile.close()
        hPairInfoFile.close()
        hOUT.close()
        return


    F1 = line1.rstrip('\n').split('\t')
    ID1 = F1[6]
    ID1 = re.sub(r'/\d$', '', ID1)
    tempID1 = ID1

    tempID2 = ""
    endFlag = 0
    for line2 in hPairInfoFile:
        F2 = line2.rstrip('\n').split('\t')
        ID2 = F2[3]
        ID2 = re.sub(r'/\d$', '', ID2)

        # treatment for the case of duplicated fastq ID (2015/11/13)
        if ID2 == tempID2: continue
        tempID2 = ID2

        while ID1 != ID2:

            line1 = hOriginalFile.readline()
            if line1  == '':
                endFlag = 1
                break

            F1 = line1.rstrip('\n').split('\t')
            ID1 = F1[6]
            ID1 = re.sub(r'/\d$', '', ID1)
            if tempID1 != "" and ID1 == tempID1: continue
            tempID1 = ID1

            if ID1 != ID2:
                print >> sys.stderr, "No pair information for %s" % ID1

        if endFlag == 1: break
 
        print >> hOUT, '\t'.join(F1[0:12]) + '\t' + F2[6] + '\t' + F2[5] + '\t' + F1[13] + '\t' + F1[14]

    hOriginalFile.close()
    hPairInfoFile.close()
    hOUT.close()



def clusterJunction(inputFilePath, outputFilePath, check_margin_size, maximum_unique_pairs):

    """
        script for merging and summarizing junction read pairs
    """
 
    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')


    mergedBedpeInfo = {}
    mergedJunction = {}
    temp_chr = None
    skip_pos = 0
    for line in hIN:

        F = line.rstrip('\n').split('\t')

        if F[0] != temp_chr: temp_chr, skip_pos = F[0], 0
        if int(F[1]) < skip_pos: continue

        match = 0
        delList = []
        for key in sorted(mergedBedpeInfo):

            tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2, inseqSize = key.split('\t')
            tids, tinseqs, tmqs1, talns1, tmqs2, talns2, tpinds, tcinds = mergedBedpeInfo[key].split('\t') 

            # the investigated key is sufficiently far from the current line in the input file and no additional line to merge is expected. therefore flush the key and information
            if F[0] != tchr1 or int(F[1]) > int(tend1) + check_margin_size:

                # obtain the most frequent junction
                junc_counter = collections.Counter(mergedJunction[key].split(';'))
                best_junc = junc_counter.most_common(1)[0][0]
                btchr1, btend1, btdir1, btchr2, btend2, btdir2, btinseq = best_junc.split(',') 
                btstart1 = str(int(btend1) - 1)
                btstart2 = str(int(btend2) - 1)


                print >> hOUT, '\t'.join([btchr1, btstart1, btend1, btchr2, btstart2, btend2, \
                                 tids, btinseq, btdir1, btdir2, tmqs1, talns1, \
                                 tmqs2, talns2, tpinds, tcinds]) + '\t' +  \
                      mergedJunction[key]

                # add to the deletion list (later the key will removed from the dictionaries)
                delList.append(key)
                continue

            else:

                # check whether the investigated key and the current line should be merged or not 
                if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2:

                    flag = 0
                    # detailed check on the junction position considering inserted sequences
                    if F[8] == "+":
                        expectedDiffSize = (int(F[2]) - int(tend1)) + (len(F[7]) - int(inseqSize))
                        if (F[9] == "+" and int(F[5]) == int(tend2) - int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(tend2) + int(expectedDiffSize)):
                            flag = 1
                    else:
                        expectedDiffSize = (int(F[2]) - int(tend1)) + (int(inseqSize) - len(F[7]))
                        if (F[9] == "+" and int(F[5]) == int(tend2) + int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(tend2) - int(expectedDiffSize)):
                            flag = 1

                    # if the junction position and direciton match
                    if flag == 1:
                
                        match = 1
                        newIds = tids + ';' + F[6]
                        newInseqs = tinseqs + ';' + F[7]
                        newMqs1 = tmqs1 + ';' + F[10]
                        newAlns1 = talns1 + ';' + F[11]
                        newMqs2 = tmqs2 + ';' + F[12]
                        newAlns2 = talns2 + ';' + F[13]
                        newPinds = tpinds + ';' + F[14]
                        newCinds = tcinds + ';' + F[15]


                        mergedBedpeInfo[key] = '\t'.join([newIds, newInseqs, newMqs1, newAlns1, newMqs2, newAlns2, newPinds, newCinds])

                        # check whether the inserted sequence should be reverse-complemented 
                        tinseq = F[7]
                        if F[7] != "---" and F[8] == F[9] and F[15] == "2":
                            # tinseq = str(Bio.Seq.Seq(F[7]).reverse_complement())
                            tinseq = utils.reverseComplement(F[7])

                        mergedJunction[key] = mergedJunction[key] + ";" + ','.join([F[0], F[2], F[8], F[3], F[5], F[9], tinseq])

        for item in delList:
            del mergedBedpeInfo[item]
            del mergedJunction[item]

        # if the current line in the input file does not match any of the pooled keys
        if match == 0:
            newKey = '\t'.join([F[0], F[1], F[2], F[3], F[4], F[5], F[8], F[9], str(len(F[7]))])
            mergedBedpeInfo[newKey] = F[6] + '\t' + F[7] + '\t' + '\t'.join(F[10:16])

            # check whether the inserted sequence should be reverse-complemented
            tinseq = F[7] 
            if F[7] != "---" and F[8] == F[9] and F[15] == "2":
                # tinseq = str(Bio.Seq.Seq(F[7]).reverse_complement())
                tinseq = utils.reverseComplement(F[7])

            mergedJunction[newKey] = ','.join([F[0], F[2], F[8], F[3], F[5], F[9], tinseq])

        if len(mergedJunction) > maximum_unique_pairs:
            print >> sys.stderr, "Exceeded maximum number of unique junction pairs at %s:%s-%s" % (F[0], F[1], F[2])
            print >> sys.stderr, "Skipp %s:%s-%s" % (F[0], F[1], str(int(F[2]) + check_margin_size))
            mergedJunction = {}
            mergedBedpeInfo = {}
            skip_pos = int(F[1]) + check_margin_size

    hIN.close()

    # last treatment
    for key in sorted(mergedBedpeInfo):

        tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2, inseqSize = key.split('\t')
        tids, tinseqs, tmqs1, talns1, tmqs2, talns2, tpinds, tcinds = mergedBedpeInfo[key].split('\t')

        # obtain the most frequent junction
        junc_counter = collections.Counter(mergedJunction[key].split(';'))
        best_junc = junc_counter.most_common(1)[0][0]
        btchr1, btend1, btdir1, btchr2, btend2, btdir2, btinseq = best_junc.split(',')
        btstart1 = str(int(btend1) - 1)
        btstart2 = str(int(btend2) - 1)

        print >> hOUT, '\t'.join([btchr1, btstart1, btend1, btchr2, btstart2, btend2, \
                         tids, btinseq, btdir1, btdir2, tmqs1, talns1, \
                         tmqs2, talns2, tpinds, tcinds]) + '\t' +  \
              mergedJunction[key]

    hOUT.close()



def parseImproperFromBam(inputBam, outputFilePath, abnormal_insert_size, min_mapping_qual, soft_clip_thres):

    """
        script for parsing improper read pairs.
        Here, we just extract on of the single-reads in the improper pairs (not read pairs),
        because the end position of the paired read cannot be extracted.
        Instead, after gathering the all the reads of improper pairs, we organize the information of improper read pairs.

        Now, the script assumes the pysam version 0.7.5. 
        This is a bit old and will be necessary to modify for the newer version of pysam 
    """

    bamfile = pysam.Samfile(inputBam, "rb")
    hOUT = open(outputFilePath, "w")

    for read in bamfile.fetch():


        # get the flag (it seems that the information about 'supplementary alignment' cannot be obtained by pysam alingmentRead class
        flags = format(int(read.flag), '#014b')[:1:-1]

        # skip unless both the reads of the pair are mapped
        if (flags[2] == "1" or flags[3] == "1"): continue

        # skip unless premaliry read
        if (flags[8] == "1" or flags[11] == "1"): continue

        # skip if the read is duplicate
        # if (flags[10] == "1"): continue

        # skip if below the minimum mapping quality
        if (read.mapq < min_mapping_qual): continue

        # skip if there is soft clipped bases in the opposite direction of alingment (in which no break point should exist). 
        if read.is_reverse is False and (read.cigar[0][0] == 4 or read.cigar[0][0] == 5) and read.cigar[0][1] >= soft_clip_thres: continue
        if read.is_reverse is True and (read.cigar[-1][0] == 4 or read.cigar[-1][0] == 5) and read.cigar[-1][1] >= soft_clip_thres: continue

        # if paired reads are not properly aligned
        # Here, not proper means
        # (I'm not sure whether the sam flags are perfect for detecting the following patterns...)
        # 1. paired reads are mapped on different chromosomes
        # 2. paired reads are mapped on the same chromosome, but very distantly ("distant" is specified by the parameter $abnormal_insert_size).
        # 3. paired reads are mapped on the same chromosome, 
        # but breaks the rule that the left-aligned read is aligned in the forward direction and the right-aligned read is in the reverse direciton.
        abnormal_pair = 0
        if read.tid != read.rnext or abs(read.isize) > abnormal_insert_size: abnormal_pair = 1

        # check for the above pattern 3
        if flags[1] == "0" and abnormal_pair == 0:
            if read.pnext - read.pos > 0 and read.is_reverse == 0 and read.mate_is_reverse == 1: continue
            if read.pos - read.pnext > 0 and read.is_reverse == 1 and read.mate_is_reverse == 0: continue
            abnormal_pair = 1


        if abnormal_pair == 1:
            seqname = (read.qname + "/1" if read.is_read1 else read.qname + "/2")  
            direction = ("-" if read.is_reverse else "+")
            print >> hOUT, seqname + '\t' + bamfile.getrname(read.tid) + '\t' + str(read.pos + 1) + '\t' + str(read.aend) + '\t' + direction + '\t' + str(read.mapq)


    bamfile.close()
    hOUT.close()



def makeImproperBedpe(inputFilePath, outputFilePath, junction_dist_margin, clipping_margin):

    ####################
    # sort according to the read ID for the later processing
    hOUT = open(outputFilePath + ".tmp1", "w")
    subprocess.call(["sort", "-k1", inputFilePath], stdout = hOUT)
    hOUT.close()
    ####################


    ####################
    # convert each improper read pair to bedpe records (with margins)
    hIN = open(outputFilePath + ".tmp1", "r")
    hOUT = open(outputFilePath + ".tmp2", "w")

    tempID, tempPairNum, tempChr, tempStart, tempEnd, tempDir, tempMapQ = "", "", "", "", "", "", ""

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        pairNum = F[0][-1:]
        F[0] = F[0][0:-2]

        if F[0] == tempID and tempPairNum == "1" and pairNum == "2":
            
            chr1, dir1, start1, end1, mapQ1, align1 = tempChr, tempDir, 0, 0, tempMapQ, tempChr + ":" + str(tempStart) + "-" + str(tempEnd) 
            if dir1 == "+":
                start1 = tempEnd - clipping_margin
                end1 = tempEnd + junction_dist_margin
            else:
                start1 = tempStart - junction_dist_margin
                end1 = tempStart + clipping_margin

            chr2, dir2, start2, end2, mapQ2, align2 = F[1], F[4], 0, 0, F[5], F[1] + ":" + F[2] + "-" + F[3]
            if dir2 == "+":
                start2 = int(F[3]) - clipping_margin
                end2 = int(F[3]) + junction_dist_margin
            else:
                start2 = int(F[2]) - junction_dist_margin
                end2 = int(F[2]) + clipping_margin
     
            if chr1 < chr2:
                print >> hOUT, '\t'.join([chr1, str(start1), str(end1), chr2, str(start2), str(end2), tempID, mapQ1 + "," + mapQ2, dir1, dir2, align1 + "," + align2])
            elif chr1 > chr2:
                print >> hOUT, '\t'.join([chr2, str(start2), str(end2), chr1, str(start1), str(end1), tempID, mapQ2 + "," + mapQ1, dir2, dir1, align2 + "," + align1])
            else:
                if start1 <= start2:
                    print >> hOUT, '\t'.join([chr1, str(start1), str(end1), chr2, str(start2), str(end2), tempID, mapQ1 + "," + mapQ2, dir1, dir2, align1 + "," + align2])
                else:
                    print >> hOUT, '\t'.join([chr2, str(start2), str(end2), chr1, str(start1), str(end1), tempID, mapQ2 + "," + mapQ1, dir2, dir1, align2 + "," + align1])


        tempID, tempPairNum, tempChr, tempStart, tempEnd, tempDir, tempMapQ = F[0], pairNum, F[1], int(F[2]), int(F[3]), F[4], F[5]

    hIN.close()
    hOUT.close()
    ####################


    ####################
    # sort according to the chromosome coordinates
    utils.sortBedpe(outputFilePath + ".tmp2", outputFilePath)
    ####################


    ####################
    # delete intermediate file
    subprocess.call(["rm", outputFilePath + '.tmp1'])
    subprocess.call(["rm", outputFilePath + '.tmp2'])
    ####################



def clusterImproperBedpe(inputFilePath, outputFilePath, check_margin_size):
   
    ####################
    # cluster and summarize improper read pair bed file
    hIN = open(inputFilePath, "r")
    hOUT = open(outputFilePath, "w")

    mergedBedpe = {}
    for line in hIN:

        F = line.rstrip('\n').split('\t')

        match = 0
        delList = []
        for key in sorted(mergedBedpe):

            tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2 = key.split('\t')
            tids, tmqs, talns = mergedBedpe[key].split('\t')

            if F[0] != tchr1 or int(F[1]) > int(tend1) + check_margin_size:
                talns_a = talns.split(';')
                talns_a_uniq = list(set(talns_a))

                if len(talns_a_uniq) >= 1:
                    
                    print >> hOUT, '\t'.join([tchr1, tstart1, tend1, tchr2, tstart2, tend2, \
                                     tids, tmqs, tdir1, tdir2, talns])
                    delList.append(key)
                    continue

            else:
                
                if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2:
                    if int(F[2]) > int(tstart1) and int(F[1]) <= int(tend1) and int(F[5]) > int(tstart2) and int(F[4]) <= int(tend2):

                        match = 1
                        newStart1 = str(max(int(tstart1), int(F[1])))  
                        newEnd1 = str(min(int(tend1), int(F[2])))
                        newStart2 = str(max(int(tstart2), int(F[4])))
                        newEnd2 = str(min(int(tend2), int(F[5])))

                        newKey = '\t'.join([tchr1, newStart1, newEnd1, tchr2, newStart2, newEnd2, tdir1, tdir2])
                        newIds = tids + ';' + F[6]
                        newMqs = tmqs + ';' + F[7]
                        newAlns = talns + ';' + F[10]
        
                        if newKey != key:    
                            delList.append(key)

                        mergedBedpe[newKey] = newIds + '\t' + newMqs + '\t' + newAlns
                        break

        for item in delList:
            del mergedBedpe[item]

        if match == 0:
            newKey = '\t'.join([F[0], F[1], F[2], F[3], F[4], F[5], F[8], F[9]])
            mergedBedpe[newKey] = F[6] + '\t' + F[7] + '\t' + F[10]

 
    for key in sorted(mergedBedpe):

        tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2 = key.split('\t')
        tids, tmqs, talns = mergedBedpe[key].split('\t')

        talns_a = talns.split(';')
        talns_a_uniq = list(set(talns_a))

        if len(talns_a_uniq) >= 1:

            print >> hOUT, '\t'.join([tchr1, tstart1, tend1, tchr2, tstart2, tend2, \
                             tids, tmqs, tdir1, tdir2, talns])


    hIN.close()
    hOUT.close()
    ####################
 

    
