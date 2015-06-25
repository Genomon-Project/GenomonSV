#! /usr/local/bin/python

import argparse
import yaml

try:
    sample_conf = yaml.load(open("sample1.yaml", "r"))
except yaml.YAMLError, exc:
    print "Error in configuration file:", exc
