#! /usr/local/bin/python

"""

    This script utilizes the SA tags (SA:Z:rname, pos, strand, CIGAR, mapQ, number of mismatch).
    The strand of the supplementary alignment information in the SA tag is determined by the orignal sequence (before taking complement).
    Therefore, please not that when the primary alignment is in the reverse direction, the sequence shown in the bam file does not match
    to the SA tags..


"""

import sys, pysam, re

inputBAM = sys.argv[1]
region = sys.argv[2]

junction_dist = 800
abnormal_insert_size = 1000
min_major_clip_size = 20
max_minor_clip_size = 15 
bamfile = pysam.Samfile(inputBAM, "rb")

regionMatch = re.search(r'(\w+):(\d+)\-(\d+)', region)
chr_region = regionMatch.group(1)
start_region = regionMatch.group(2)
end_region = regionMatch.group(3)

SAre = re.compile('(\w+),(\d+),([\-\+]),(\w+),(\d+)')
cigarMDRe = re.compile('(\d+)([MD])')
cigarHIMSRe = re.compile('(\d+)([HIMS])')
cigarHSRe_right = re.compile('(\d+)([HS])$')
cigarHSRe_left = re.compile('^(\d+)([HS])')

for read in bamfile.fetch(chr_region, int(start_region), int(end_region)):

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
                print '\t'.join([juncChr_current, str(juncPos_current - 1), str(juncPos_current), juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), \
                                 read.qname + ("/1" if flags[6] == "1" else "/2"), juncSurplus, juncDir_current, juncDir_SA, \
                                 str(read.mapq), coverRegion_current + "," + coverRegion_SA, chr_pair + ":" + str(pos_pair), str(juncType), "1"])
            else: 
                print '\t'.join([juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), juncChr_current, str(juncPos_current - 1), str(juncPos_current), \
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
                print '\t'.join([juncChr_current, str(juncPos_current - 1), str(juncPos_current), juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), \
                                 read.qname + ("/1" if flags[6] == "1" else "/2"), juncSurplus, juncDir_current, juncDir_SA, \
                                 str(read.mapq), coverRegion_current + "," + coverRegion_SA, chr_pair + ":" + str(pos_pair), str(juncType), "1"])
            else:                
                print '\t'.join([juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), juncChr_current, str(juncPos_current - 1), str(juncPos_current), \
                                 read.qname + ("/1" if flags[6] == "1" else "/2"), juncSurplus, juncDir_SA, juncDir_current, \
                                 str(read.mapq), coverRegion_current + "," + coverRegion_SA, chr_pair + ":" + str(pos_pair), str(juncType), "2"])

