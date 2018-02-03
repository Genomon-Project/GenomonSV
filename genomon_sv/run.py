#!/usr/bin/env python


import sys, argparse, subprocess, os, multiprocessing
# import config 
import utils
import parseFunction
import filterFunction
import mergeFunction
import annotationFunction

def genomonSV_parse(args):

    """
        Genomon SV: parsing breakpoint containing read pairs and improperly aligned read pairs
    """

    ####################
    # check the existence of input ba
    if not os.path.exists(args.bam_file):
        raise ValueError('No file: ' + args.bam_file)

    ####################
    # make output directories
    utils.make_directory(os.path.dirname(args.output_prefix))
    ####################

    
    ####################
    # parse breakpoint containing read pairs from input bam files
    utils.processingMessage("parsing breakpoint containing read pairs from bam file") 
    parseFunction.parseJunctionFromBam(args.bam_file, args.output_prefix + ".junction.unsort.txt", args.junction_min_mapping_qual,
                                       args.junction_abnormal_insert_size, args.junction_min_major_clipping_size, args.junction_max_minor_clipping_size)

    utils.processingMessage("sorting parsed breakpoint containing read pairs")
    utils.sortBedpe(args.output_prefix + ".junction.unsort.txt", args.output_prefix + ".junction.sort.txt")

    utils.processingMessage("getting start positions of paired-end reads")
    parseFunction.getPairStartPos(args.output_prefix + ".junction.sort.txt", args.output_prefix + ".junction.pairStart.bed")

    utils.processingMessage("compressing start position infomation file")
    utils.compress_index_bed(args.output_prefix + ".junction.pairStart.bed", args.output_prefix + ".junction.pairStart.bed.gz")


    utils.processingMessage("getting covered regions of paired-end reads from bam file")
    parseFunction.getPairCoverRegionFromBam(args.bam_file, args.output_prefix + ".junction.pairCoverage.txt", args.output_prefix + ".junction.pairStart.bed.gz")

    utils.processingMessage("adding information of covered regions of paired-end reads")
    parseFunction.addPairCoverRegionFromBam(args.output_prefix + ".junction.sort.txt", args.output_prefix + ".junction.sort.withPair.txt", args.output_prefix + ".junction.pairCoverage.txt")

    utils.processingMessage("clustering breakpoint containing read pairs")
    parseFunction.clusterJunction(args.output_prefix + ".junction.sort.withPair.txt", args.output_prefix + ".junction.clustered.bedpe.unsort", \
                                  args.junction_check_margin_size, args.junction_check_maximum_unique_pairs)

    utils.processingMessage("sorting clustered breakpoint containing read pairs")
    utils.sortBedpe(args.output_prefix + ".junction.clustered.bedpe.unsort", args.output_prefix + ".junction.clustered.bedpe")

    utils.processingMessage("compressing clustered breakpoint containing read pairs")
    utils.compress_index_bed(args.output_prefix + ".junction.clustered.bedpe", args.output_prefix + ".junction.clustered.bedpe.gz")

    if args.debug == False:
        subprocess.call(["rm", args.output_prefix + ".junction.unsort.txt"])
        subprocess.call(["rm", args.output_prefix + ".junction.sort.txt"])
        subprocess.call(["rm", args.output_prefix + ".junction.pairStart.bed.gz"])
        subprocess.call(["rm", args.output_prefix + ".junction.pairStart.bed.gz.tbi"])
        subprocess.call(["rm", args.output_prefix + ".junction.sort.withPair.txt"])
        subprocess.call(["rm", args.output_prefix + ".junction.pairCoverage.txt"])
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.bedpe.unsort"])
        subprocess.call(["rm", args.output_prefix + ".junction.pairStart.bed"])
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.bedpe"])
 
    ####################
    # improper read pairs

    # parse potentially improper read pairs from input bam files
    utils.processingMessage("parsing improperly aligned read pairs from bam file")
    parseFunction.parseImproperFromBam(args.bam_file, args.output_prefix + ".improper.unsort.txt",
                                       args.improper_abnormal_insert_size, args.improper_min_mapping_qual, args.improper_max_clipping_size)

    # create and organize bedpe file integrating pair information
    utils.processingMessage("sorting improperly aligned read pairs") 
    parseFunction.makeImproperBedpe(args.output_prefix + ".improper.unsort.txt", args.output_prefix + ".improper.bedpe", 
                                    args.junction_dist_margin, args.junction_opposite_dist_margin)

    # cluster read pairs possibly representing the same junction
    utils.processingMessage("clustering improperly aligned read pairs")
    parseFunction.clusterImproperBedpe(args.output_prefix + ".improper.bedpe", args.output_prefix + ".improper.clustered.unsort.bedpe", args.improper_check_margin_size)

    utils.processingMessage("sorting clustered improperly aligned read pairs")
    utils.sortBedpe(args.output_prefix + ".improper.clustered.unsort.bedpe", args.output_prefix + ".improper.clustered.bedpe")

    utils.processingMessage("compressing clustered improperly aligned read pairs")
    utils.compress_index_bed(args.output_prefix + ".improper.clustered.bedpe", args.output_prefix + ".improper.clustered.bedpe.gz")

    if args.debug == False:
        subprocess.call(["rm", args.output_prefix + ".improper.unsort.txt"])
        subprocess.call(["rm", args.output_prefix + ".improper.bedpe"])
        subprocess.call(["rm", args.output_prefix + ".improper.clustered.unsort.bedpe"])
        subprocess.call(["rm", args.output_prefix + ".improper.clustered.bedpe"])

    ####################

def genomonSV_filt(args):

    ####################
    # file existence check
    if not os.path.exists(args.bam_file):
        raise ValueError('No file: ' + args.bam_file)

    if not os.path.exists(args.output_prefix + ".junction.clustered.bedpe.gz"):
        raise ValueError('No file: ' + args.output_prefix + ".junction.clustered.bedpe.gz")

    if not os.path.exists(args.output_prefix + ".junction.clustered.bedpe.gz.tbi"):
        raise ValueError('No file: ' + args.output_prefix + ".junction.clustered.bedpe.gz.tbi")

    if not os.path.exists(args.output_prefix + ".improper.clustered.bedpe.gz"):
        raise ValueError('No file: ' + args.output_prefix + ".improper.clustered.bedpe.gz")

    if not os.path.exists(args.output_prefix + ".improper.clustered.bedpe.gz.tbi"):
        raise ValueError('No file: ' + args.output_prefix + ".improper.clustered.bedpe.gz.tbi")

    if args.matched_control_bam != "" and not os.path.exists(args.matched_control_bam):
        raise ValueError('No file: ' + args.matched_control_bam)

    if args.non_matched_control_junction != "" and not os.path.exists(args.non_matched_control_junction):
        raise ValueError('No file: ' + args.non_matched_control_junction)

    if not os.path.exists(args.reference_genome):
        raise ValueError('No file: ' + args.reference_genome)

    # if not os.path.exists(args.annotation_dir + "/gene.bed.gz"):
    #     raise ValueError('No file: ' + args.annotation_dir + "/gene.bed.gz")

    # if not os.path.exists(args.annotation_dir + "/exon.bed.gz"):
    #     raise ValueError('No file: ' + args.annotation_dir + "/exon.bed.gz")

    if args.thread_num == 1:
        filterFunction.genomon_sv_filt_main(args.output_prefix, args)
    else:
        thread_num_mod = filterFunction.partition_junction(args.output_prefix, args.thread_num)

        jobs = []
        for i in range(1, thread_num_mod + 1):
            proc = multiprocessing.Process(target = filterFunction.genomon_sv_filt_main, \
                                           args = (args.output_prefix + ".thread_" + str(i), args, " (thread " + str(i) + ')'))
            jobs.append(proc)
            proc.start()

        for i in range(0, thread_num_mod):
            jobs[i].join()

        header_flag = 0
        with open(args.output_prefix + ".genomonSV.result.txt", 'w') as hout:
            for i in range(1, thread_num_mod + 1):
                with open(args.output_prefix + ".thread_" + str(i) + ".genomonSV.result.txt", 'r') as hin:
                    header = hin.readline().rstrip('\n')
                    if header_flag == 0: 
                        print >> hout, header
                        header_flag = 1
                    for line in hin:
                       print >> hout, line.rstrip('\n') 
 
        for i in range(1, thread_num_mod + 1):
            subprocess.check_call(["rm", "-rf", args.output_prefix + ".thread_" + str(i) + ".junction.clustered.bedpe.gz"])
            subprocess.check_call(["rm", "-rf", args.output_prefix + ".thread_" + str(i) + ".junction.clustered.bedpe.gz.tbi"])
            subprocess.check_call(["rm", "-rf", args.output_prefix + ".thread_" + str(i) + ".genomonSV.result.txt"])


        
    """ 
    ####################
    utils.processingMessage("filtering by # of breakpoint containing read pairs and variant sizes")
    filterFunction.filterJuncNumAndSize(args.output_prefix + ".junction.clustered.bedpe.gz",
                                        args.output_prefix + ".junction.clustered.filt1.bedpe",
                                        args.min_junc_num, args.min_sv_size, args.min_inversion_size)

    utils.processingMessage("filtering by nonmatched control panel")
    filterFunction.filterNonMatchControl(args.output_prefix + ".junction.clustered.filt1.bedpe",
                                         args.output_prefix + ".junction.clustered.filt2.bedpe",
                                         args.non_matched_control_junction,
                                         args.matched_control_label,
                                         args.control_panel_num_thres, args.control_panel_check_margin)

    utils.processingMessage("incorporating improperly alinged read pair infomation")
    filterFunction.addImproperInfo(args.output_prefix + ".junction.clustered.filt2.bedpe",
                                   args.output_prefix + ".junction.clustered.filt3.bedpe",
                                   args.output_prefix + ".improper.clustered.bedpe.gz")

    utils.processingMessage("filtering by sizes of covered regions, mapping quality and # of support read pairs")
    filterFunction.filterMergedJunc(args.output_prefix + ".junction.clustered.filt3.bedpe",
                                    args.output_prefix + ".junction.clustered.filt4.bedpe",
                                    args.min_support_num, args.min_mapping_qual, args.min_overhang_size)

    utils.processingMessage("filtering too close candidates")
    filterFunction.removeClose(args.output_prefix + ".junction.clustered.filt4.bedpe",
                               args.output_prefix + ".junction.clustered.filt5.bedpe",
                               args.close_check_margin, args.close_check_thres)

    utils.processingMessage("performing realignments")
    filterFunction.validateByRealignment(args.output_prefix + ".junction.clustered.filt5.bedpe",
                    args.output_prefix + ".junction.clustered.filt6.bedpe",
                    args.bam_file, args.matched_control_bam, args.reference_genome, args.blat_option,
                    args.short_tandem_reapeat_thres, args.max_depth, args.search_length, args.search_margin, 
                    args.split_refernece_thres, args.validate_sequence_length)

    utils.processingMessage("filtering allele frequencies, Fisher's exact test p-values and # of support read pairs")
    filterFunction.filterNumAFFis(args.output_prefix + ".junction.clustered.filt6.bedpe", 
                                  args.output_prefix + ".junction.clustered.filt7.bedpe",
                                  args.matched_control_bam,
                                  args.min_tumor_variant_read_pair, args.min_tumor_allele_freq, 
                                  args.max_control_variant_read_pair, args.max_control_allele_freq,
                                  args.max_fisher_pvalue)

    utils.processingMessage("adding annotation")
    annotationFunction.addAnnotation(args.output_prefix + ".junction.clustered.filt7.bedpe",
                                     args.output_prefix + ".genomonSV.result.txt",
                                     args.genome_id, args.grc)

    if args.debug == False:
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.filt1.bedpe"])
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.filt2.bedpe"])
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.filt3.bedpe"])
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.filt4.bedpe"])
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.filt5.bedpe"])
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.filt6.bedpe"])
        subprocess.call(["rm", args.output_prefix + ".junction.clustered.filt7.bedpe"])

    ####################
    """


def genomonSV_merge(args):

    """
    script for merging clustered junction data for creating nonmatched normal control data
    the first input the path for the individual junction list file to be merged.
    """

    with open(args.control_info_file, 'r') as hin:
        for line in hin:
            label, output_prefix = line.rstrip('\n').split('\t')
            if not os.path.exists(output_prefix + ".junction.clustered.bedpe.gz"):
                raise ValueError('No file: ' + output_prefix + ".junction.clustered.bedpe.gz")                


    utils.make_directory(os.path.dirname(args.merge_output_file))
    
    if os.path.exists(args.merge_output_file + ".temp"):
        print >> sys.stderr, "Remove existing intermediate file " + args.merge_output_file + ".temp" 
        os.remove(args.merge_output_file + ".temp")

    with open(args.control_info_file, 'r') as hin:
        for line in hin:
            label, output_prefix = line.rstrip('\n').split('\t')
            utils.processingMessage("extracting information of " + label)
            mergeFunction.simplifyJunc(output_prefix + ".junction.clustered.bedpe.gz", args.merge_output_file + ".temp", label)

    
    utils.processingMessage("sorting the aggregated junction file")
    utils.sortBedpe(args.merge_output_file + ".temp", args.merge_output_file + ".temp.sort")

    utils.processingMessage("merging the same junction in the aggregated junction file")
    mergeFunction.organizeControl(args.merge_output_file + ".temp.sort", args.merge_output_file + ".temp.merged", args.merge_check_margin_size)

    utils.processingMessage("sorting the merged junction file")
    utils.sortBedpe(args.merge_output_file + ".temp.merged", args.merge_output_file + ".temp.merged.sort")

    utils.processingMessage("compressing the merged junction file")
    utils.compress_index_bed(args.merge_output_file + ".temp.merged.sort", args.merge_output_file)


    if args.debug == False:
        subprocess.call(["rm", args.merge_output_file + ".temp"])
        subprocess.call(["rm", args.merge_output_file + ".temp.sort"])
        subprocess.call(["rm", args.merge_output_file + ".temp.merged"])
        subprocess.call(["rm", args.merge_output_file + ".temp.merged.sort"])

