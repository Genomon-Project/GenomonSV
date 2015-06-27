#!/usr/bin/env python

"""

    Genomon SV filt: extracting candidates of structural variations from clustered breakpoint-containing and improperly aligned read pairs

"""

import sys
import argparse
import config 
import utils
import parseFunction
import filterFunction 

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
    # file existence check

 
    ####################
    outputPrefix = sampleConf["output_dir"] + "/" + sampleConf["label"]

    filterFunction.filterJuncNumAndSize(outputPrefix + ".junction.clustered.bedpe.gz",
                                        outputPrefix + ".junction.clustered.filt1.bedpe",
                                        paramConf["filterCondition"])

    ####################


if __name__ == "__main__":
    genomonSV_parse()
