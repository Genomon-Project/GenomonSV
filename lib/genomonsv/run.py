#!/usr/bin/env python


import sys, argparse, subprocess, os
import config 
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
    utils.make_directory(os.path.dirname(args.parse_output_prefix))
    ####################

    
    ####################
    # parse breakpoint containing read pairs from input bam files
    utils.processingMessage("parsing breakpoint containing read pairs from bam file") 
    parseFunction.parseJunctionFromBam(args.bam_file, args.parse_output_prefix + ".junction.unsort.txt",
                                       args.junction_abnormal_insert_size, args.junction_min_major_clipping_size, args.junction_max_minor_clipping_size)

    utils.processingMessage("sorting parsed breakpoint containing read pairs")
    utils.sortBedpe(args.parse_output_prefix + ".junction.unsort.txt", args.parse_output_prefix + ".junction.sort.txt")

    utils.processingMessage("getting start positions of paired-end reads")
    parseFunction.getPairStartPos(args.parse_output_prefix + ".junction.sort.txt", args.parse_output_prefix + ".junction.pairStart.bed")

    utils.processingMessage("compressing start position infomation file")
    utils.compress_index_bed(args.parse_output_prefix + ".junction.pairStart.bed", args.parse_output_prefix + ".junction.pairStart.bed.gz")


    utils.processingMessage("getting covered regions of paired-end reads from bam file")
    parseFunction.getPairCoverRegionFromBam(args.bam_file, args.parse_output_prefix + ".junction.pairCoverage.txt", args.parse_output_prefix + ".junction.pairStart.bed.gz")

    utils.processingMessage("adding information of covered regions of paired-end reads")
    parseFunction.addPairCoverRegionFromBam(args.parse_output_prefix + ".junction.sort.txt", args.parse_output_prefix + ".junction.sort.withPair.txt", args.parse_output_prefix + ".junction.pairCoverage.txt")

    utils.processingMessage("clustering breakpoint containing read pairs")
    parseFunction.clusterJunction(args.parse_output_prefix + ".junction.sort.withPair.txt", args.parse_output_prefix + ".junction.clustered.bedpe.unsort", args.junction_check_margin_size)

    utils.processingMessage("sorting clustered breakpoint containing read pairs")
    utils.sortBedpe(args.parse_output_prefix + ".junction.clustered.bedpe.unsort", args.parse_output_prefix + ".junction.clustered.bedpe")

    utils.processingMessage("compressing clustered breakpoint containing read pairs")
    utils.compress_index_bed(args.parse_output_prefix + ".junction.clustered.bedpe", args.parse_output_prefix + ".junction.clustered.bedpe.gz")

    if args.debug == False:
        subprocess.call(["rm", args.parse_output_prefix + ".junction.unsort.txt"])
        subprocess.call(["rm", args.parse_output_prefix + ".junction.sort.txt"])
        subprocess.call(["rm", args.parse_output_prefix + ".junction.pairStart.bed.gz"])
        subprocess.call(["rm", args.parse_output_prefix + ".junction.pairStart.bed.gz.tbi"])
        subprocess.call(["rm", args.parse_output_prefix + ".junction.sort.withPair.txt"])
        subprocess.call(["rm", args.parse_output_prefix + ".junction.pairCoverage.txt"])
        subprocess.call(["rm", args.parse_output_prefix + ".junction.clustered.bedpe.unsort"])
        subprocess.call(["rm", args.parse_output_prefix + ".junction.pairStart.bed"])
        subprocess.call(["rm", args.parse_output_prefix + ".junction.clustered.bedpe"])
 
    ####################
    # improper read pairs

    # parse potentially improper read pairs from input bam files
    utils.processingMessage("parsing improperly aligned read pairs from bam file")
    parseFunction.parseImproperFromBam(args.bam_file, args.parse_output_prefix + ".improper.unsort.txt",
                                       args.improper_abnormal_insert_size, args.improper_min_mapping_qual, args.improper_max_clipping_size)

    # create and organize bedpe file integrating pair information
    utils.processingMessage("sorting improperly aligned read pairs") 
    parseFunction.makeImproperBedpe(args.parse_output_prefix + ".improper.unsort.txt", args.parse_output_prefix + ".improper.bedpe", 
                                    args.junction_dist_margin, args.junction_opposite_dist_margin_margin)

    # cluster read pairs possibly representing the same junction
    utils.processingMessage("clustering improperly aligned read pairs")
    parseFunction.clusterImproperBedpe(args.parse_output_prefix + ".improper.bedpe", args.parse_output_prefix + ".improper.clustered.unsort.bedpe", args.improper_check_margin_size)

    utils.processingMessage("sorting clustered improperly aligned read pairs")
    utils.sortBedpe(args.parse_output_prefix + ".improper.clustered.unsort.bedpe", args.parse_output_prefix + ".improper.clustered.bedpe")

    utils.processingMessage("compressing clustered improperly aligned read pairs")
    utils.compress_index_bed(args.parse_output_prefix + ".improper.clustered.bedpe", args.parse_output_prefix + ".improper.clustered.bedpe.gz")

    if args.debug == False:
        subprocess.call(["rm", args.parse_output_prefix + ".improper.unsort.txt"])
        subprocess.call(["rm", args.parse_output_prefix + ".improper.bedpe"])
        subprocess.call(["rm", args.parse_output_prefix + ".improper.clustered.unsort.bedpe"])
        subprocess.call(["rm", args.parse_output_prefix + ".improper.clustered.bedpe"])

    ####################

def genomonSV_filt(args):


    ####################
    # load config files
    global sampleConf
    sampleConf = config.sample_yaml_config_parse(args.sampleInfoFile, "filt")

    global paramConf
    paramConf = config.param_yaml_contig_parse(args.paramInfoFile, "filt")
    ####################


    ####################
    # file existence check

 
    ####################
    outputPrefix = sampleConf["target"]["path_to_output_dir"] + "/" + sampleConf["target"]["label"]
    matchedControlFlag = sampleConf["matched_control"]["use"]
    matchedControlBam = sampleConf["matched_control"]["path_to_bam"] if matchedControlFlag == True else None
    
    use_non_matched_control_panel = sampleConf["non_matched_control_panel"]["use"]
    data_path_non_matched_control_panel = sampleConf["non_matched_control_panel"]["data_path"] if use_non_matched_control_panel == True else "---"
    matched_control_label_non_matched_control_panel = sampleConf["non_matched_control_panel"]["matched_control_label"] if use_non_matched_control_panel == True else "---"

    utils.processingMessage("filtering by # of breakpoint containing read pairs and variant sizes")
    filterFunction.filterJuncNumAndSize(outputPrefix + ".junction.clustered.bedpe.gz",
                                        outputPrefix + ".junction.clustered.filt1.bedpe",
                                        paramConf["filter_condition"])


    utils.processingMessage("filtering by nonmatched control panel")
    filterFunction.filterNonMatchControl(outputPrefix + ".junction.clustered.filt1.bedpe",
                                         outputPrefix + ".junction.clustered.filt2.bedpe",
                                         use_non_matched_control_panel,
                                         data_path_non_matched_control_panel,
                                         matched_control_label_non_matched_control_panel,
                                         paramConf["filter_condition"])
    # subprocess.call(["cp", outputPrefix + ".junction.clustered.filt1.bedpe", outputPrefix + ".junction.clustered.filt2.bedpe"])

    utils.processingMessage("incorporating improperly alinged read pair infomation")
    filterFunction.addImproperInfo(outputPrefix + ".junction.clustered.filt2.bedpe",
                                   outputPrefix + ".junction.clustered.filt3.bedpe",
                                   outputPrefix + ".improper.clustered.bedpe.gz")

    utils.processingMessage("filtering by sizes of covered regions, mapping quality and # of support read pairs")
    filterFunction.filterMergedJunc(outputPrefix + ".junction.clustered.filt3.bedpe",
                                    outputPrefix + ".junction.clustered.filt4.bedpe",
                                    paramConf["filter_condition"])

    utils.processingMessage("filtering too close candidates")
    filterFunction.removeClose(outputPrefix + ".junction.clustered.filt4.bedpe",
                               outputPrefix + ".junction.clustered.filt5.bedpe",
                               paramConf["filter_condition"])

    utils.processingMessage("performing realignments")
    filterFunction.validateByRealignment(outputPrefix + ".junction.clustered.filt5.bedpe",
                    outputPrefix + ".junction.clustered.filt6.bedpe",
                    sampleConf["target"]["path_to_bam"],
                    matchedControlBam,
                    paramConf["software"]["blat"] + " " + paramConf["software"]["blat_option"],
                    matchedControlFlag,
                    paramConf["realignment_validation_condition"])

    utils.processingMessage("filtering allele frequencies, Fisher's exact test p-values and # of support read pairs")
    filterFunction.filterNumAFFis(outputPrefix + ".junction.clustered.filt6.bedpe", 
                                  outputPrefix + ".junction.clustered.filt7.bedpe",
                                  matchedControlFlag,
                                  paramConf["realignment_validation_condition"])

    utils.processingMessage("adding annotation")
    annotationFunction.addAnnotation(outputPrefix + ".junction.clustered.filt7.bedpe",
                                     outputPrefix + ".genomonSV.result.txt",
                                     paramConf["annotation"])

    if paramConf["debug_mode"] == False:
        subprocess.call(["rm", outputPrefix + ".junction.clustered.filt1.bedpe"])
        subprocess.call(["rm", outputPrefix + ".junction.clustered.filt2.bedpe"])
        subprocess.call(["rm", outputPrefix + ".junction.clustered.filt3.bedpe"])
        subprocess.call(["rm", outputPrefix + ".junction.clustered.filt4.bedpe"])
        subprocess.call(["rm", outputPrefix + ".junction.clustered.filt5.bedpe"])
        subprocess.call(["rm", outputPrefix + ".junction.clustered.filt6.bedpe"])
        subprocess.call(["rm", outputPrefix + ".junction.clustered.filt7.bedpe"])

    ####################


def genomonSV_merge(args):

    """
    script for merging clustered junction data for creating nonmatched normal control data
    the first input the path for the individual junction list file to be merged.
    """

    with open(args.control_info_file, 'r') as hin:
        for line in hin:
            label, parse_output_prefix = line.rstrip('\n').split('\t')
            if not os.path.exists(parse_output_prefix + ".junction.clustered.bedpe.gz"):
                raise ValueError('No file: ' + parse_output_prefix + ".junction.clustered.bedpe.gz")                


    utils.make_directory(os.path.dirname(args.merge_output_file))
    
    if os.path.exists(args.merge_output_file + ".temp"):
        print >> sys.stderr, "Remove existing intermediate file " + args.merge_output_file + ".temp" 
        os.remove(args.merge_output_file + ".temp")

    with open(args.control_info_file, 'r') as hin:
        for line in hin:
            label, parse_output_prefix = line.rstrip('\n').split('\t')
            utils.processingMessage("extracting information of " + label)
            mergeFunction.simplifyJunc(parse_output_prefix + ".junction.clustered.bedpe.gz", args.merge_output_file + ".temp", label)

    
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

