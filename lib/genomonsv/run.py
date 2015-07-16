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
    # load config files
    global sampleConf
    sampleConf = config.sample_yaml_config_parse(args.sampleInfoFile, "parse")

    global paramConf
    paramConf = config.param_yaml_contig_parse(args.paramInfoFile, "parse")
    ####################


    ####################
    # make output directories
    utils.make_directory(sampleConf["target"]["path_to_output_dir"])
    ####################

    
    ####################
    outputPrefix = sampleConf["target"]["path_to_output_dir"] + "/" + sampleConf["target"]["label"]

    # parse breakpoint containing read pairs from input bam files
    utils.processingMessage("parsing breakpoint containing read pairs from bam file") 
    parseFunction.parseJunctionFromBam(sampleConf["target"]["path_to_bam"], 
                                       outputPrefix + ".junction.unsort.txt", 
                                       paramConf["parse_junction_condition"])

    utils.processingMessage("sorting parsed breakpoint containing read pairs")
    utils.sortBedpe(outputPrefix + ".junction.unsort.txt",
                    outputPrefix + ".junction.sort.txt")

    utils.processingMessage("getting start positions of paired-end reads")
    parseFunction.getPairStartPos(outputPrefix + ".junction.sort.txt",
                                  outputPrefix + ".junction.pairStart.bed")

    utils.processingMessage("compressing start position infomation file")
    utils.compress_index_bed(outputPrefix + ".junction.pairStart.bed",
                             outputPrefix + ".junction.pairStart.bed.gz",
                             paramConf["software"]["bgzip"], paramConf["software"]["tabix"])


    utils.processingMessage("getting covered regions of paired-end reads from bam file")
    parseFunction.getPairCoverRegionFromBam(sampleConf["target"]["path_to_bam"], 
                                            outputPrefix + ".junction.pairCoverage.txt",
                                            outputPrefix + ".junction.pairStart.bed.gz")

    utils.processingMessage("adding information of covered regions of paired-end reads")
    parseFunction.addPairCoverRegionFromBam(outputPrefix + ".junction.sort.txt",
                                            outputPrefix + ".junction.sort.withPair.txt",
                                            outputPrefix + ".junction.pairCoverage.txt")

    utils.processingMessage("clustering breakpoint containing read pairs")
    parseFunction.clusterJunction(outputPrefix + ".junction.sort.withPair.txt", 
                                  outputPrefix + ".junction.clustered.bedpe.unsort",
                                  paramConf["cluster_junction_condition"])

    utils.processingMessage("sorting clustered breakpoint containing read pairs")
    utils.sortBedpe(outputPrefix + ".junction.clustered.bedpe.unsort", outputPrefix + ".junction.clustered.bedpe")

    utils.processingMessage("compressing clustered breakpoint containing read pairs")
    utils.compress_index_bed(outputPrefix + ".junction.clustered.bedpe",
                             outputPrefix + ".junction.clustered.bedpe.gz",
                             paramConf["software"]["bgzip"], paramConf["software"]["tabix"])

    if paramConf["debug_mode"] == False:
        subprocess.call(["rm", outputPrefix + ".junction.unsort.txt"])
        subprocess.call(["rm", outputPrefix + ".junction.sort.txt"])
        subprocess.call(["rm", outputPrefix + ".junction.pairStart.bed.gz"])
        subprocess.call(["rm", outputPrefix + ".junction.pairStart.bed.gz.tbi"])
        subprocess.call(["rm", outputPrefix + ".junction.sort.withPair.txt"])
        subprocess.call(["rm", outputPrefix + ".junction.pairCoverage.txt"])
        subprocess.call(["rm", outputPrefix + ".junction.clustered.bedpe.unsort"])
        subprocess.call(["rm", outputPrefix + ".junction.pairStart.bed"])
        subprocess.call(["rm", outputPrefix + ".junction.clustered.bedpe"])
 
    ####################
    # improper read pairs

    # parse potentially improper read pairs from input bam files
    utils.processingMessage("parsing improperly aligned read pairs from bam file")
    parseFunction.parseImproperFromBam(sampleConf["target"]["path_to_bam"],
                         outputPrefix + ".improper.unsort.txt",
                         paramConf["parse_improper_condition"])

    # create and organize bedpe file integrating pair information
    utils.processingMessage("sorting improperly aligned read pairs") 
    parseFunction.makeImproperBedpe(outputPrefix + ".improper.unsort.txt",
                                    outputPrefix + ".improper.bedpe",
                                    paramConf["cluster_improper_condition"])

    # cluster read pairs possibly representing the same junction
    utils.processingMessage("clustering improperly aligned read pairs")
    parseFunction.clusterImproperBedpe(outputPrefix + ".improper.bedpe",
                                       outputPrefix + ".improper.clustered.unsort.bedpe",
                                       paramConf["cluster_improper_condition"])

    utils.processingMessage("sorting clustered improperly aligned read pairs")
    utils.sortBedpe(outputPrefix + ".improper.clustered.unsort.bedpe",
                    outputPrefix + ".improper.clustered.bedpe")

    utils.processingMessage("compressing clustered improperly aligned read pairs")
    utils.compress_index_bed(outputPrefix + ".improper.clustered.bedpe",
                             outputPrefix + ".improper.clustered.bedpe.gz",
                             paramConf["software"]["bgzip"], paramConf["software"]["tabix"])

    if paramConf["debug_mode"] == False:
        subprocess.call(["rm", outputPrefix + ".improper.unsort.txt"])
        subprocess.call(["rm", outputPrefix + ".improper.bedpe"])
        subprocess.call(["rm", outputPrefix + ".improper.clustered.unsort.bedpe"])
        subprocess.call(["rm", outputPrefix + ".improper.clustered.bedpe"])

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


    utils.processingMessage("filtering by # of breakpoint containing read pairs and variant sizes")
    filterFunction.filterJuncNumAndSize(outputPrefix + ".junction.clustered.bedpe.gz",
                                        outputPrefix + ".junction.clustered.filt1.bedpe",
                                        paramConf["filter_condition"])


    if sampleConf["non_matched_control_panel"]["use"] == True:
        utils.processingMessage("filtering by nonmatched control panel")
        filterFunction.filterNonMatchControl(outputPrefix + ".junction.clustered.filt1.bedpe",
                                             outputPrefix + ".junction.clustered.filt2.bedpe",
                                             sampleConf["non_matched_control_panel"]["data_path"],
                                             sampleConf["non_matched_control_panel"]["matched_control_label"],
                                             paramConf["filter_condition"])
    else:
        subprocess.call(["cp", outputPrefix + ".junction.clustered.filt1.bedpe", outputPrefix + ".junction.clustered.filt2.bedpe"])

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
                    sampleConf["matched_control"]["path_to_bam"],
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

    ####################
    # load config files
    controlConf = config.control_yaml_config_parse(args.controlPathInfoFile)

    outputFilePath = args.outputFilePath

    paramConf = config.param_yaml_contig_parse(args.paramInfoFile, "merge")
    ####################

    
    if os.path.exists(outputFilePath + ".temp"):
        print >> sys.stderr, "Remove existing intermediate file " + outputFilePath + ".temp"
        os.remove(outputFilePath + ".temp")

    for label in sorted(controlConf):
        utils.processingMessage("extracting information of " + controlConf[label])
        mergeFunction.simplifyJunc(controlConf[label], outputFilePath + ".temp", label)

    utils.processingMessage("sorting the aggregated junction file")
    utils.sortBedpe(outputFilePath + ".temp",
                    outputFilePath + ".temp.sort")

    utils.processingMessage("merging the same junction in the aggregated junction file")
    mergeFunction.organizeControl(outputFilePath + ".temp.sort",
                                  outputFilePath + ".temp.merged",
                                  paramConf["control_merge_condition"])

    utils.processingMessage("sorting the merged junction file")
    utils.sortBedpe(outputFilePath + ".temp.merged",
                    outputFilePath + ".temp.merged.sort")

    utils.processingMessage("compressing the merged junction file")
    utils.compress_index_bed(outputFilePath + ".temp.merged.sort",
                             outputFilePath,
                             paramConf["software"]["bgzip"], paramConf["software"]["tabix"])

    if paramConf["debug_mode"] == False:
        subprocess.call(["rm", outputFilePath + ".temp"])
        subprocess.call(["rm", outputFilePath + ".temp.sort"])
        subprocess.call(["rm", outputFilePath + ".temp.merged"])


