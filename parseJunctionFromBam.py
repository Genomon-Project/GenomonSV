#! /usr/local/bin/python

import sys, pysam, re

inputBAM = sys.argv[1]
region = sys.argv[2]

junction_dist = 800
abnormal_insert_size = 1000
min_clip_size_for_primary = 20

bamfile = pysam.Samfile(inputBAM, "rb")

regionMatch = re.search(r'(\w+):(\d+)\-(\d+)', region)
chr_region = regionMatch.group(1)
start_region = regionMatch.group(2)
end_region = regionMatch.group(3)

SAre = re.compile('(\w+),(\d+),([\-\+]),(\w+),(\d+)')
cigarMDRe = re.compile('(\d+)([MD])')
cigarHIMSRe = re.compile('(\d+)([HIMS])')
cigarHSRe_right = re.compile('(\d+)([HIMS])$')
cigarHSRe_left = re.compile('^(\d+)([HIMS])')

for read in bamfile.fetch(chr_region, int(start_region), int(end_region)):

    if read.qname == "ST-E00129:110:H0E3PALXX:1:1111:23451:14863":
        print read.qname

    # get the flag information
    flags = format(int(read.flag), "#014b")[:1:-1]

    # skip if either of the read pair is unmapped
    if flags[2] == "1" or flags[3] == "1": continue

    # skip supplementary alignment
    if flags[8] == "1" or flags[11] == "1": continue

    # skip if the read aligned to hs37d5"
    # (in the future, this step will be replaced to some more sophisticated way;
    # (e.g., the user can input the target chromosomes and ignore if the read is aligned to non-target chromosomes, and so on..
    if read.tid == "hs37d5" or read.rnext == "hs37d5": continue

    # no clipping
    if len(read.cigar) == 1: continue

    # get the clipping size in the both side
    left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
    right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)

    if left_clipping < min_clip_size_for_primary and right_clipping < min_clip_size_for_primary: continue

    # for comparing with the previous script (this will removed soon)
    if left_clipping > 0 and right_clipping > 0: continue

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

    if right_clipping > min_clip_size_for_primary:

        chr_SA, pos_SA, dir_SA, cigar_SA = SA_str.group(1), int(SA_str.group(2)), SA_str.group(3), SA_str.group(4)
        clipLen_current = right_clipping
        alignmentSize_current = read.alen
        readLength_current = read.rlen

        juncChr_current = chr_current
        juncPos_current = pos_current + alignmentSize_current - 1
        juncDir_current = "+"
        juncPos_SA = chr_SA
        coverRegion_current = chr_current + ":" + str(pos_current) + "-" + str(juncPos_current)
        juncChr_SA = chr_SA

        expected_clipLen_SA = readLength_current - clipLen_current
        expected_clipDir_SA = ("-" if dir_current == dir_SA else "+")

        alignmentSize_SA = 0
        for item in cigarMDRe.finditer(cigar_SA):
            alignmentSize_SA += int(item.group(1))        

        # get the soft clipping information on the supplementary alignment
        right_clipping_SA = 0 
        tempMatch = cigarHSRe_right.match(cigar_SA)
        if tempMatch is not None: right_clipping_SA = int(tempMatch.group(1))

        left_clipping_SA = 0 
        tempMatch = cigarHSRe_left.match(cigar_SA)
        if tempMatch is not None: left_clipping_SA = int(tempMatch.group(1))


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
                validFlag = 1

            if dir_SA == "-" and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                clipLen_SA = left_clipping_SA 
                juncDir_SA = "-"
                juncPos_SA = pos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                coverRegion_SA = chr_SA + ":" + str(juncPos_SA) + "-" + str(juncPos_SA + alignmentSize_SA - 1)
                validFlag = 1


        if dir_current == "+" and chr_SA == chr_pair:

            if dir_SA == "-" and dir_pair == "+" and 0 <= pos_SA - pos_pair < abnormal_insert_size and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                clipLen_SA = right_clipping_SA
                juncDir_SA = "+"
                juncPos_SA = pos_SA + alignmentSize_SA - 1
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                coverRegion_SA = chr_SA + ":" + str(pos_SA) + "-" + str(juncPos_SA) 
                validFlag = 1

            if dir_SA == "+" and dir_pair == "-" and 0 <= pos_pair - pos_SA < abnormal_insert_size and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                clipLen_SA = left_clipping_SA
                juncDir_SA = "-"
                juncPos_SA = pos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                coverRegion_SA = chr_SA + ":" + str(juncPos_SA) + "-" + str(juncPos_SA + alignmentSize_SA - 1)
                validFlag = 1


        if validFlag == 1:

            juncSurplus = "---"
            if clipLen_SA > expected_clipLen_SA and readLength_current == len(read.seq):
                surPlus_start = read.alen - clipLen_current
                surPlus_end = surPlus_start + clipLen_SA - expected_clipLen_SA
                juncSurplus = read.seq[surPlus_start:surPlus_end]

            print '\t'.join([juncChr_current, str(juncPos_current - 1), str(juncPos_current), juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), \
                             read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq), juncDir_current, juncDir_SA, \
                             coverRegion_current + "," + coverRegion_SA, juncSurplus, chr_pair + ":" + str(pos_pair)])

 

    if left_clipping > min_clip_size_for_primary:

        chr_SA, pos_SA, dir_SA, cigar_SA = SA_str.group(1), int(SA_str.group(2)), SA_str.group(3), SA_str.group(4)
        clipLen_current = left_clipping
        alignmentSize_current = read.alen
        readLength_current = read.rlen
 
        juncChr_current = chr_current
        juncPos_current = pos_current
        juncDir_current = "-"
        juncPos_SA = chr_SA
        coverRegion_current = chr_current + ":" + str(pos_current) + "-" + str(pos_current + alignmentSize_current - 1)
        juncChr_SA = chr_SA

        expected_clipLen_SA = readLength_current - clipLen_current
        expected_clipDir_SA = ("+" if dir_current == dir_SA else "-")

        alignmentSize_SA = 0
        for item in cigarMDRe.finditer(cigar_SA):
            alignmentSize_SA += int(item.group(1))

        # get the soft clipping information on the supplementary alignment
        right_clipping_SA = 0
        tempMatch = cigarHSRe_right.match(cigar_SA)
        if tempMatch is not None: right_clipping_SA = int(tempMatch.group(1))

        left_clipping_SA = 0
        tempMatch = cigarHSRe_left.match(cigar_SA)
        if tempMatch is not None: left_clipping_SA = int(tempMatch.group(1))


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
                validFlag = 1

            if dir_SA == "-" and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                clipLen_SA = left_clipping_SA
                juncDir_SA = "-"
                juncPos_SA = pos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                coverRegion_SA = chr_SA + ":" + str(juncPos_SA) + "-" + str(juncPos_SA + alignmentSize_SA - 1)
                validFlag = 1


        if dir_current == "-" and chr_SA == chr_pair:

           if dir_SA == "-" and dir_pair == "+" and 0 <= pos_SA - pos_pair < abnormal_insert_size and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                clipLen_SA = right_clipping_SA
                juncDir_SA = "+"
                juncPos_SA = pos_SA + alignmentSize_SA - 1
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                coverRegion_SA = chr_SA + ":" + str(pos_SA) + "-" + str(juncPos_SA)
                validFlag = 1

           if dir_SA == "+" and dir_pair == "-" and 0 <= pos_pair - pos_SA < abnormal_insert_size and expected_clipDir_SA == "-" and left_clipping_SA:
                clipLen_SA = left_clipping_SA
                juncDir_SA = "-"
                juncPos_SA = pos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                coverRegion_SA = chr_SA + ":" + str(juncPos_SA) + "-" + str(juncPos_SA + alignmentSize_SA - 1)
                validFlag = 1


        if validFlag == 1:

            juncSurplus = "---"
            if clipLen_SA > expected_clipLen_SA and readLength_current == len(read.seq):
                surPlus_start = read.alen - clipLen_current
                surPlus_end = surPlus_start + clipLen_SA - expected_clipLen_SA
                juncSurplus = read.seq[surPlus_start:surPlus_end]

            print '\t'.join([juncChr_current, str(juncPos_current - 1), str(juncPos_current), juncChr_SA, str(juncPos_SA - 1), str(juncPos_SA), \
                             read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq), juncDir_current, juncDir_SA, \
                             coverRegion_current + "," + coverRegion_SA, juncSurplus, chr_pair + ":" + str(pos_pair)])


