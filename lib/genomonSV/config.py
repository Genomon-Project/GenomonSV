#!/usr/bin/env python

import yaml

def sample_yaml_config_parse(filePath):

    try:
        with open(filePath, 'r') as fIN:
            sampleConf = yaml.load(fIN)
    except yaml.YAMLError, exc:
        print "Error in sample information file:", exc

    # check whether necessary information is provided in right way
    # check the existence of input files
    #
        
    return sampleConf


def param_yaml_contig_parse(filePath):

    try:
        with open(filePath, 'r') as fIN:
            paramConf = yaml.load(fIN)
    except yaml.YAMLError, exc:
        print "Error in sample information file:", exc

    # check whether necessary information is provided in right way
    #
    #
    return paramConf

