#!/usr/bin/env python


import sys
import argparse
import config 
import utils
import parseFunction
import filterFunction

def genomonSV_parse(args):

    """
        Genomon SV: parsing breakpoint containing read pairs and improperly aligned read pairs
    """

    ####################
    # load config files
    global sampleConf
    sampleConf = config.sample_yaml_config_parse(args.sampleInfoFile)

    global paramConf
    paramConf = config.param_yaml_contig_parse(args.paramInfoFile)
    ####################


    ####################
    # make output directories
    utils.make_directory(sampleConf["target"]["outputDir"])
    ####################

    
    ####################
    outputPrefix = sampleConf["target"]["outputDir"] + "/" + sampleConf["target"]["label"]

    # parse breakpoint containing read pairs from input bam files
    parseFunction.parseJunctionFromBam(sampleConf["target"]["path_to_bam"], 
                                       outputPrefix + ".junction.unsort.txt", 
                                       paramConf["parseJunctionCondition"])

    utils.sortBedpe(outputPrefix + ".junction.unsort.txt",
                    outputPrefix + ".junction.sort.txt")

    parseFunction.getPairStartPos(outputPrefix + ".junction.sort.txt",
                                  outputPrefix + ".junction.pairStart.bed")

    utils.compress_index_bed(outputPrefix + ".junction.pairStart.bed",
                             outputPrefix + ".junction.pairStart.bed.gz",
                             paramConf["software"]["bgzip"], paramConf["software"]["tabix"])



    parseFunction.getPairCoverRegionFromBam(sampleConf["target"]["path_to_bam"], 
                                            outputPrefix + ".junction.pairCoverage.txt",
                                            outputPrefix + ".junction.pairStart.bed.gz")


    parseFunction.addPairCoverRegionFromBam(outputPrefix + ".junction.sort.txt",
                                            outputPrefix + ".junction.sort.withPair.txt",
                                            outputPrefix + ".junction.pairCoverage.txt")

    parseFunction.clusterJunction(outputPrefix + ".junction.sort.withPair.txt", 
                                  outputPrefix + ".junction.clustered.bedpe.unsort",
                                  paramConf["clusterJunctionCondition"])

    utils.sortBedpe(outputPrefix + ".junction.clustered.bedpe.unsort", outputPrefix + ".junction.clustered.bedpe")

    utils.compress_index_bed(outputPrefix + ".junction.clustered.bedpe",
                             outputPrefix + ".junction.clustered.bedpe.gz",
                             paramConf["software"]["bgzip"], paramConf["software"]["tabix"])

    ####################
    # improper read pairs

    # parse potentially improper read pairs from input bam files
    parseFunction.parseImproperFromBam(sampleConf["target"]["path_to_bam"],
                         outputPrefix + ".improper.unsort.txt",
                         paramConf["parseImproperCondition"])

    # create and organize bedpe file integrating pair information 
    parseFunction.makeImproperBedpe(outputPrefix + ".improper.unsort.txt",
                                    outputPrefix + ".improper.bedpe",
                                    paramConf["clusterImproperCondition"])

    # cluster read pairs possibly representing the same junction
    parseFunction.clusterImproperBedpe(outputPrefix + ".improper.bedpe",
                                       outputPrefix + ".improper.clustered.bedpe",
                                       paramConf["clusterImproperCondition"])

    utils.compress_index_bed(outputPrefix + ".improper.clustered.bedpe",
                             outputPrefix + ".improper.clustered.bedpe.gz",
                             paramConf["software"]["bgzip"], paramConf["software"]["tabix"])
    ####################


def genomonSV_filt(args):

    """
        Genomon SV filt: extracting candidates of structural variations from clustered breakpoint-containing and improperly aligned read pairs
    """

    """
    ####################
    # parse arguments
    parser = argparse.ArgumentParser(description = "Parse bam file for breakpoint and improper read pairs")
    parser.add_argument("sampleInfoFile", metavar = "sample.yaml", type = str,
                        help = "input sample information file (yaml format)")
    parser.add_argument("paramInfoFile", metavar = "param.yaml", type = str,
                        help = "parameter information file (yaml format)")
    args = parser.parse_args()
    ####################
    """

    ####################
    # load config files
    global sampleConf
    sampleConf = config.sample_yaml_config_parse(args.sampleInfoFile)

    global paramConf
    paramConf = config.param_yaml_contig_parse(args.paramInfoFile)
    ####################


    ####################
    # file existence check

 
    ####################
    outputPrefix = sampleConf["target"]["outputDir"] + "/" + sampleConf["target"]["label"]

    filterFunction.filterJuncNumAndSize(outputPrefix + ".junction.clustered.bedpe.gz",
                                        outputPrefix + ".junction.clustered.filt1.bedpe",
                                        paramConf["filterCondition"])

    filterFunction.filterNonMatchControl(outputPrefix + ".junction.clustered.filt1.bedpe",
                                         outputPrefix + ".junction.clustered.filt2.bedpe",
                                         controlFile, matchedNormal, Params)
    ####################


