#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import genomon_sv 
from check_download import *

class TestFilt(unittest.TestCase):

    def setUp(self):
        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")
 
        self.parser = genomon_sv.arg_parser.create_parser()


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        shutil.copyfile(cur_dir + "/data/parse/5929_tumor.junction.clustered.bedpe.gz", tmp_dir + "/5929_tumor.junction.clustered.bedpe.gz")
        shutil.copyfile(cur_dir + "/data/parse/5929_tumor.junction.clustered.bedpe.gz.tbi", tmp_dir + "/5929_tumor.junction.clustered.bedpe.gz.tbi")
        shutil.copyfile(cur_dir + "/data/parse/5929_tumor.improper.clustered.bedpe.gz", tmp_dir + "/5929_tumor.improper.clustered.bedpe.gz")
        shutil.copyfile(cur_dir + "/data/parse/5929_tumor.improper.clustered.bedpe.gz.tbi", tmp_dir + "/5929_tumor.improper.clustered.bedpe.gz.tbi")

        tumor_bam = cur_dir + "/data/bam/5929_tumor.markdup.bam"
        output_prefix = tmp_dir + "/5929_tumor"
        control_bam = cur_dir + "/data/bam/5929_control.markdup.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
 
        output_file = tmp_dir + "/5929_tumor.genomonSV.result.txt"
        answer_file = cur_dir + "/data/parse/5929_tumor.genomonSV.result.txt"

        print ' '.join(["filt", tumor_bam, output_prefix, ref_genome, "--grc", "--matched_control_bam", control_bam])
        args = self.parser.parse_args(["filt", tumor_bam, output_prefix, ref_genome, "--grc", "--matched_control_bam", control_bam])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


    def test2(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        shutil.copyfile(cur_dir + "/data/parse/5929_tumor.junction.clustered.bedpe.gz", tmp_dir + "/5929_tumor.junction.clustered.bedpe.gz")
        shutil.copyfile(cur_dir + "/data/parse/5929_tumor.junction.clustered.bedpe.gz.tbi", tmp_dir + "/5929_tumor.junction.clustered.bedpe.gz.tbi")
        shutil.copyfile(cur_dir + "/data/parse/5929_tumor.improper.clustered.bedpe.gz", tmp_dir + "/5929_tumor.improper.clustered.bedpe.gz")
        shutil.copyfile(cur_dir + "/data/parse/5929_tumor.improper.clustered.bedpe.gz.tbi", tmp_dir + "/5929_tumor.improper.clustered.bedpe.gz.tbi")

        tumor_bam = cur_dir + "/data/bam/5929_tumor.markdup.bam"
        output_prefix = tmp_dir + "/5929_tumor"
        control_bam = cur_dir + "/data/bam/5929_control.markdup.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"

        output_file = tmp_dir + "/5929_tumor.genomonSV.result.txt"
        answer_file = cur_dir + "/data/parse/5929_tumor.genomonSV.result.txt"

        print ' '.join(["filt", tumor_bam, output_prefix, ref_genome, "--grc", "--matched_control_bam", control_bam, "--thread_num", "4"])
        args = self.parser.parse_args(["filt", tumor_bam, output_prefix, ref_genome, "--grc", "--matched_control_bam", control_bam, "--thread_num", "4"])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        # shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()

