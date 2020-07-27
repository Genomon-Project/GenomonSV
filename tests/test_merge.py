#! /usr/bin/env python

from __future__ import print_function
import unittest
import os, glob, tempfile, shutil, filecmp
import genomon_sv 
import gzip

class TestMerge(unittest.TestCase):

    def setUp(self):
        self.parser = genomon_sv.arg_parser.create_parser()

    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        all_junc_file = glob.glob(cur_dir + "/data/parse/*.junction.clustered.bedpe.gz")
        with open(tmp_dir + "/control_info.txt", 'w') as hout:
            for junc_file in sorted(all_junc_file):
                sample = os.path.basename(junc_file).replace(".junction.clustered.bedpe.gz", '')
                prefix = junc_file.replace(".junction.clustered.bedpe.gz", '')
                print("%s\t%s" % (sample, prefix), file = hout)

        control_list_file = tmp_dir + "/control_info.txt"
        output_file = tmp_dir + "/merge_control.bedpe.gz"
        answer_file = cur_dir + "/data/merge/merge_control.bedpe.gz"

        args = self.parser.parse_args(["merge", control_list_file, output_file])
        args.func(args)

        tmp_answer_file = tmp_dir+"/answer_merge_control.bedpe"
        with gzip.open(answer_file,"rt") as hin:
            with open(tmp_answer_file,"w") as hout:
                shutil.copyfileobj(hin,hout)

        tmp_output_file = tmp_dir+"/output_merge_control.bedpe"
        with gzip.open(output_file,"rt") as hin:
            with open(tmp_output_file,"w") as hout:
                shutil.copyfileobj(hin,hout)

        self.assertTrue(filecmp.cmp(tmp_output_file, tmp_answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()

