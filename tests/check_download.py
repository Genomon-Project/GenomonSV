#! /usr/bin/env python

import os, urllib2

def check_download(url, output_path):

    if not os.path.exists(os.path.dirname(output_path)):
       os.makedirs(os.path.dirname(output_path))

    if not os.path.exists(output_path):
        ufile = urllib2.urlopen(url)
        with open(output_path, 'w') as hout:
            for x in ufile:
                hout.write(x)

