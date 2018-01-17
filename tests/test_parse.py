#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import genomon_sv 

class TestParse(unittest.TestCase):

    def setUp(self):
        self.parser = genomon_sv.arg_parser.create_parser()

    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/bam/5929_tumor.markdup.bam"

        output_prefix = tmp_dir + "/5929_tumor"
        output_file1 = tmp_dir + "/5929_tumor.junction.clustered.bedpe.gz"
        output_file2 = tmp_dir + "/5929_tumor.junction.clustered.bedpe.gz.tbi"
        output_file3 = tmp_dir + "/5929_tumor.improper.clustered.bedpe.gz"
        output_file4 = tmp_dir + "/5929_tumor.improper.clustered.bedpe.gz.tbi"

        answer_file1 = cur_dir + "/data/parse/5929_tumor.junction.clustered.bedpe.gz"
        answer_file2 = cur_dir + "/data/parse/5929_tumor.junction.clustered.bedpe.gz.tbi"
        answer_file3 = cur_dir + "/data/parse/5929_tumor.improper.clustered.bedpe.gz"
        answer_file4 = cur_dir + "/data/parse/5929_tumor.improper.clustered.bedpe.gz.tbi"
 
        args = self.parser.parse_args(["parse", input_file, output_prefix])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file1, answer_file1, shallow=False))
        self.assertTrue(filecmp.cmp(output_file2, answer_file2, shallow=False))
        self.assertTrue(filecmp.cmp(output_file3, answer_file3, shallow=False))
        self.assertTrue(filecmp.cmp(output_file4, answer_file4, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()

