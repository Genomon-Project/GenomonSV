#!/usr/bin/env python

import sys, pysam, subprocess
import annot_utils.gene, annot_utils.exon

def addAnnotation(inputFilePath, outputFilePath, genome_id, is_grc):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    annot_utils.gene.make_gene_info(outputFilePath + ".tmp.refGene.bed.gz", "refseq", genome_id, is_grc, False)
    annot_utils.exon.make_exon_info(outputFilePath + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, False)

    gene_bed = outputFilePath + ".tmp.refGene.bed.gz" 
    exon_bed = outputFilePath + ".tmp.refExon.bed.gz"

    gene_tb = pysam.TabixFile(gene_bed)
    exon_tb = pysam.TabixFile(exon_bed)

    print >> hOUT, '\t'.join(["Chr_1", "Pos_1", "Dir_1", "Chr_2", "Pos_2", "Dir_2", "Inserted_Seq", "Variant_Type", \
                             "Gene_1", "Gene_2", "Exon_1", "Exon_2", "Num_Tumor_Ref_Read_Pair", "Num_Tumor_Var_Read_Pair", \
                             "Tumor_VAF", "Num_Control_Ref_Read_Pair", "Num_Control_Var_Read_Pair", "Control_VAF", \
                             "Minus_Log_Fisher_P_value", "Non-Matched_Control_Sample_With_Max_Junction", "Num_Max_Non-Matched_Control_Junction", \
                             "Max_Over_Hang_1", "Max_Over_Hang_2"])
 
    for line in hIN:

        F = line.rstrip('\n').split('\t')

        # tumorAF = 0
        # if float(F[7]) + float(F[8]) > 0: tumorAF = float(F[8]) / (float(F[7]) + float(F[8]))
        # normalAF = 0
        # if float(F[9]) + float(F[10]) > 0: normalAF = float(F[10]) / (float(F[9]) + float(F[10]))

        ##########
        # check gene annotation for the side 1  
        tabixErrorFlag = 0
        try:
            records = gene_tb.fetch(F[0], int(F[1]) - 1, int(F[1]))
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1

        gene1 = [];
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                gene1.append(record[3])

        if len(gene1) == 0: gene1.append("---")
        gene1 = list(set(gene1))
        ##########

        ##########
        # check gene annotation for the side 2
        tabixErrorFlag = 0
        try:
            records = gene_tb.fetch(F[3], int(F[4]) - 1, int(F[4]))
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1

        gene2 = [];
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                gene2.append(record[3])

        if len(gene2) == 0: gene2.append("---")
        gene2 = list(set(gene2))
        ##########

        ##########
        # check exon annotation for the side 1
        tabixErrorFlag = 0
        try:
            records = exon_tb.fetch(F[0], int(F[1]) - 1, int(F[1]))
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1

        exon1 = [];
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                exon1.append(record[3])

        if len(exon1) == 0: exon1.append("---")
        exon1 = list(set(exon1))
        ##########

        ##########
        # check exon annotation for the side 2
        tabixErrorFlag = 0
        try:
            records = exon_tb.fetch(F[3], int(F[4]) - 1, int(F[4]))
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1
       
        exon2 = [];
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
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

        print >> hOUT, '\t'.join(F[0:7]) + '\t' + SVtype + '\t' + ';'.join(gene1) + '\t' + ';'.join(gene2) + '\t' + ';'.join(exon1) + '\t' + ';'.join(exon2) + '\t' + \
                       '\t'.join(F[7:]) 
#               '\t'.join(F[7:11]) + '\t' + str(round(tumorAF, 4)) + '\t' + str(round(normalAF, 4)) + '\t' + str(round(float(F[11]), 4))

    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.refGene.bed.gz"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.refGene.bed.gz.tbi"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.refExon.bed.gz"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.refExon.bed.gz.tbi"])


    hIN.close()
    hOUT.close()
    gene_tb.close()
    exon_tb.close()


