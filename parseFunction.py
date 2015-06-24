#! /usr/local/bin/python

"""
    functions for parsing breakpoint containing read pairs and improperly aligned read pairs
"""

import pysam, re, subprocess

def parseJunctionFromBam(inputBAM, outputFilePath, Params):

    """
    This function utilizes the SA tags (SA:Z:rname, pos, strand, CIGAR, mapQ, number of mismatch).
    The strand of the supplementary alignment information in the SA tag is determined by the orignal sequence (before taking complement).
    Therefore, please not that when the primary alignment is in the reverse direction, the sequence shown in the bam file does not match
    to the SA tags..
    """

    abnormal_insert_size = Params["abnormal_insert_size"]
    min_major_clip_size = Params["min_major_clip_size"]
    max_minor_clip_size = Params["max_minor_clip_size"]

    bamfile = pysam.Samfile(inputBAM, "rb")
    hOUT = open(outputFilePath, "w")
 
    SAre = re.compile('(\w+),(\d+),([\-\+]),(\w+),(\d+)')
    cigarMDRe = re.compile('(\d+)([MD])')
    cigarHIMSRe = re.compile('(\d+)([HIMS])')
    cigarHSRe_right = re.compile('(\d+)([HS])$')
    cigarHSRe_left = re.compile('^(\d+)([HS])')


    # maybe add the regional extraction of bam files
    for read in bamfile.fetch()):

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



def parseImproperFromBam(inputBam, outputFilePath, Params):

    """
        script for parsing improper read pairs.
        Here, we just extract on of the single-reads in the improper pairs (not read pairs),
        because the end position of the paired read cannot be extracted.
        Instead, after gathering the all the reads of improper pairs, we organize the information of improper read pairs.

        Now, the script assumes the pysam version 0.7.5. 
        This is a bit old and will be necessary to modify for the newer version of pysam 
    """

    abnormal_insert_size = Params["abnormal_insert_size"]
    min_mapping_qual = Params["min_mapping_qual"]
    soft_clip_thres = Params["soft_clip_thres"]

    bamfile = pysam.Samfile(inputBAM, "rb")
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


def sortJunction(inputFilePath, outputFilePath):

    hOUT = open(outputFile, "w")
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", inputFilePath], stdout = hOUT)
    hOUT.close()


def makePairStartBed(inputFilePath, outputFilePath, Param):

    ####################
    # obtain the position information about the pair read from the junction file
    hIN = open(inputFilePath, "r")
    hOUT = open(outputFilePath + '.tmp', "w")

    num = 1
    reChrPos = re.compile('^(\w+):(\d+)')
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

        print chr + '\t' + str(int(pos) - 1) + '\t' + pos + '\t' + ID + '\t' + str(num)

        num = num + 1

    hIN.close()
    hOUT.close() 
    ####################
 

    ####################
    # sort the temporary file
    hOUT = open(outputFilePath, "w")
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", inputFilePath], stdout = hOUT)
    hOUT.close()
    ####################


    ####################
    # compress by bgzip
    hOUT = open(outputFilePath + ".gz", "w")
    subprocess.call(["bgzip", "-f", outputFilePath], stdout = hOUT)
    hOUT.close()
    ####################


    ####################
    # index by tabix
    subprocess.call(["tabix", "-p", "bed", outputFilePath + ".gz"])
    ####################


    ####################
    # delete intermediate file
    subprocess.call(["rm", outputFilePath + '.tmp'])



# def makeImproperBedpe(inputFilePath, outputFilePath, Param):




    
