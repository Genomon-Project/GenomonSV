#! /usr/local/bin/python

"""

    Genomon SV: parsing breakpoint containing read pairs and improperly aligned read pairs

"""

import sys
import argparse
import yaml

def genomonSV_parse():

    ####################
    # parse arguments
    parser = argparse.ArgumentParser(description = "Parse bam file for breakpoint and improper read pairs")
    parser.add_argument("sampleInfoFile", metavar = "sample.yaml", type=argparse.FileType('r'), 
                        help = "input sample information file (yaml format)")
    parser.add_argument("paramInfoFile", metavar = "param.yaml", type=argparse.FileType('r'),
                        help = "parameter information file (yaml format)")
    args = parser.parse_args()
    ####################

    ####################
    # load config files
    try:
        sampleConf = yaml.load(args.sampleInfoFile)
    except yaml.YAMLError, exc:
        print "Error in sample information file:", exc

    try:
        paramConf = yaml.load(args.paramInfoFile)
    except yaml.YAMLError, exc:
        print "Error in sample information file:", exc

    # checking the config files (hopefully, create conf.py to check each key
    ####################



    # check the existence of input files

    # make output directories


    # parse breakpoint containing read pairs from input bam files
    # parseJunctionFromBam()

    # parse improper read pairs from input bam files
    # parseImproperFromBam()

    # sort and compress the parsed breakpoint containing read pairs
    # sortJunction()

    # merge and cluster parsed improper read pairs from input bam files
    # clusterImproper()

    # merge and cluster parsed breakpoint containing read pairs
    # clusterJunction()





if __name__ == "__main__":
    genomonSV_parse()
