#!/usr/bin/env python

import yaml, os

def sample_yaml_config_parse(filePath, method):

    try:
        with open(filePath, 'r') as fIN:
            sampleConf = yaml.load(fIN)
    except yaml.YAMLError, exc:
        print "Error in sample information file:", exc

    # check whether necessary information is provided in right way
    # target:
    if not "target" in sampleConf:
        raise ValueError('No key: target. Please confirm sample.yaml config file.')

    if not "label" in sampleConf["target"]:
        raise ValueError('No key: target->label. Please confirm sample.yaml config file.')
 
    if not "path_to_bam" in sampleConf["target"]:
        raise ValueError('No key: target->path_to_bam. Please confirm sample.yaml config file.')

    if not "path_to_output_dir" in sampleConf["target"]:
        raise ValueError('No key: target->path_to_output_dir. Please confirm sample.yaml config file.')

    if not os.path.exists(sampleConf["target"]["path_to_bam"]):
        raise ValueError('No file: ' + sampleConf["target"]["path_to_bam"])


    if method == "filt":

        if not "matched_control" in sampleConf:
            raise ValueError('No key: matched_control. Please confirm sample.yaml config file.')

        if not "use" in sampleConf["matched_control"]:
            raise ValueError('No key: matched_control->use. Please confirm sample.yaml config file.')

        if not isinstance(sampleConf["matched_control"]["use"], bool):
            raise ValueError('The key matched_control->use should be boolean (True or False). Please confirm sample.yaml config file.')

        if sampleConf["matched_control"]["use"] == True:
            if not "path_to_bam" in sampleConf["matched_control"]:
                raise ValueError('No key: matched_control->path_to_bam. Please confirm sample.yaml config file.')

            if not os.path.exists(sampleConf["matched_control"]["path_to_bam"]):
                raise ValueError('No file: ' + sampleConf["matched_control"]["path_to_bam"])


        if not "non_matched_control_panel" in sampleConf:
            raise ValueError('No key: non_matched_control_panel. Please confirm sample.yaml config file.')

        if not "use" in sampleConf["non_matched_control_panel"]:
            raise ValueError('No key: non_matched_control_panel->use. Please confirm sample.yaml config file.')

        if not isinstance(sampleConf["non_matched_control_panel"]["use"], bool):
            raise ValueError('The key non_matched_control_panel->use should be boolean (True or False). Please confirm sample.yaml config file.')

        if sampleConf["non_matched_control_panel"]["use"] == True:
            if not "matched_control_label" in sampleConf["non_matched_control_panel"]:
                raise ValueError('No key: non_matched_control_panel->data_path. Please confirm sample.yaml config file.')

            if not "data_path" in sampleConf["non_matched_control_panel"]:
                raise ValueError('No key: matched_control->data_path. Please confirm sample.yaml config file.')

            if not os.path.exists(sampleConf["non_matched_control_panel"]["data_path"]):
                raise ValueError('No file: ' + sampleConf["non_matched_control_panel"]["data_path"])


    return sampleConf


def param_yaml_contig_parse(filePath, method):

    try:
        with open(filePath, 'r') as fIN:
            paramConf = yaml.load(fIN)
    except yaml.YAMLError, exc:
        print "Error in sample information file:", exc

    # check whether necessary information is provided in right way
    #
    #
    return paramConf



def control_yaml_config_parse(filePath):

    try:
        with open(filePath, 'r') as fIN:
            controlConf = yaml.load(fIN)
    except yaml.YAMLError, exc:
        print "Error in sample information file:", exc

    for label in controlConf:
        
        # check the exisitence of files
        if os.path.exists(controlConf[label]) == False:
            sys.exit(controlConf[label] + " in " + filePath + "does not exists!")

    return controlConf


