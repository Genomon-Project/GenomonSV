#! /usr/local/bin/python

"""

    Genomon SV: parsing breakpoint containing read pairs and improperly aligned read pairs

"""

import sys
import argparse
import config 
import utils
import parseFunction

def genomonSV_parse():

    ####################
    # parse arguments
    parser = argparse.ArgumentParser(description = "Parse bam file for breakpoint and improper read pairs")
    parser.add_argument("sampleInfoFile", metavar = "sample.yaml", type = str,
                        help = "input sample information file (yaml format)")
    parser.add_argument("paramInfoFile", metavar = "param.yaml", type = str,
                        help = "parameter information file (yaml format)")
    args = parser.parse_args()
    ####################


    ####################
    # load config files
    global sampleConf
    sampleConf = config.sample_yaml_config_parse(args.sampleInfoFile)

    global paramConf
    paramConf = config.param_yaml_contig_parse(args.paramInfoFile)
    ####################


    ####################
    # make output directories
    utils.make_directory(sampleConf["output_dir"])
    ####################



    # parse breakpoint containing read pairs from input bam files
    parseFunction.parseJunctionFromBam(sampleConf["path_to_bam"], 
                                       sampleConf["output_dir"] + "/" + sampleConf["label"] + ".junction.unsort.txt", 
                                       paramConf["parseJunctionCondition"])

    utils.sortBedpe(sampleConf["output_dir"] + "/" + sampleConf["label"] + ".junction.unsort.txt",
                    sampleConf["output_dir"] + "/" + sampleConf["label"] + ".junction.sort.txt")

    # parse improper read pairs from input bam files
    parseFunction.parseImproperFromBam(sampleConf["path_to_bam"],
                         sampleConf["output_dir"] + "/" + sampleConf["label"] + ".improper.unsort.txt",
                         paramConf["parseImproperCondition"])

    # parseFunction.makeImproperBedpe()
    # sort and compress the parsed breakpoint containing read pairs
    # sortJunction()

    # merge and cluster parsed improper read pairs from input bam files
    # clusterImproper()

    # merge and cluster parsed breakpoint containing read pairs
    # clusterJunction()





if __name__ == "__main__":
    genomonSV_parse()
