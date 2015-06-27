#! /usr/local/bin/python

"""
    script for parsing improper read pairs.
    Here, we just extract on of the single-reads in the improper pairs (not read pairs),
    because the end position of the paired read cannot be extracted.
    Instead, after gathering the all the reads of improper pairs, we organize the information of improper read pairs.

    Now, the script assumes the pysam version 0.7.5. 
    This is a bit old and will be necessary to modify for the newer version of pysam 
"""

import sys, re, pysam

inputBAM = sys.argv[1]
region = sys.argv[2]


junction_dist = 800
abnormal_insert_size = 2000
min_mapping_qual = 30
soft_clip_thres = 5

bamfile = pysam.Samfile(inputBAM, "rb")

regionMatch = re.search(r'(\w+):(\d+)\-(\d+)', region)
chr_region = regionMatch.group(1)
start_region = regionMatch.group(2)
end_region = regionMatch.group(3)

for read in bamfile.fetch(chr_region, int(start_region), int(end_region)):


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
        print seqname + '\t' + bamfile.getrname(read.tid) + '\t' + str(read.pos + 1) + '\t' + str(read.aend) + '\t' + direction + '\t' + str(read.mapq)



