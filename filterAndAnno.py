#! /usr/local/bin/python

import sys, tabix

inputFile = sys.argv[1]
geneFile = sys.argv[2]
exonFile = sys.argv[3]

hIN = open(inputFile, 'r')
gene_tb = tabix.open(geneFile)
exon_tb = tabix.open(exonFile)

for line in hIN:
    F = line.rstrip('\n').split('\t')

    tumorAF = 0
    if float(F[7]) + float(F[8]) > 0: tumorAF = float(F[8]) / (float(F[7]) + float(F[8]))     
    normalAF = 0
    if float(F[9]) + float(F[10]) > 0: normalAF = float(F[10]) / (float(F[9]) + float(F[10]))

    if int(F[8]) < 3: continue
    if tumorAF < 0.05: continue
    if normalAF > 0.2: continue
    if float(F[11]) < 1: continue
     
    ##########
    # check gene annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = gene_tb.query(F[0], int(F[1]) - 1, int(F[1]))
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    gene1 = [];
    if tabixErrorFlag == 0:
        for record in records:
            gene1.append(record[3])

    if len(gene1) == 0: gene1.append("---")
    gene1 = list(set(gene1))
    ##########

    ##########
    # check gene annotation for the side 2
    tabixErrorFlag = 0
    try:
        records = gene_tb.query(F[3], int(F[4]) - 1, int(F[4]))
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    gene2 = [];
    if tabixErrorFlag == 0:
        for record in records:
            gene2.append(record[3])

    if len(gene2) == 0: gene2.append("---")
    gene2 = list(set(gene2))
    ##########

    ##########
    # check exon annotation for the side 1
    tabixErrorFlag = 0
    try:
        records = exon_tb.query(F[0], int(F[1]) - 1, int(F[1]))
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    exon1 = [];
    if tabixErrorFlag == 0:
        for record in records:
            exon1.append(record[3])

    if len(exon1) == 0: exon1.append("---")
    exon1 = list(set(exon1))
    ##########

    ##########
    # check exon annotation for the side 2
    tabixErrorFlag = 0
    try:
        records = exon_tb.query(F[3], int(F[4]) - 1, int(F[4]))
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1
   
    exon2 = [];
    if tabixErrorFlag == 0:
        for record in records:
            exon2.append(record[3])

    if len(exon2) == 0: exon2.append("---")
    exon2 = list(set(exon2))
    ##########

    if F[0] != F[3]:
        SVtype = "translocation"
    elif F[2] == "+" and F[5] == "-":
        SVtype = "deletion"
    elif F[2] == "-" and F[5] == "+":
        SVtype = "tandem_duplication"
    else:
        SVtype = "inversion"

    print '\t'.join(F[0:7]) + '\t' + SVtype + '\t' + ';'.join(gene1) + '\t' + ';'.join(gene2) + '\t' + ';'.join(exon1) + '\t' + ';'.join(exon2) + '\t' + \
          '\t'.join(F[7:11]) + '\t' + str(round(tumorAF, 4)) + '\t' + str(round(normalAF, 4)) + '\t' + str(round(float(F[11]), 4))
  

hIN.close()



