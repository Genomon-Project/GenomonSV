#! /usr/bin/env python

import unittest
import os, glob, tempfile, shutil, filecmp
import genomon_sv 

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
                print >> hout, "%s\t%s" % (sample, prefix)

        control_list_file = tmp_dir + "/control_info.txt"
        output_file = tmp_dir + "/merge_control.bedpe.gz"
        answer_file = cur_dir + "/data/merge/merge_control.bedpe.gz"

        args = self.parser.parse_args(["merge", control_list_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()

