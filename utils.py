import os, subprocess

def make_directory(inputDir):
    """
    make input directory if it does not exist.
    """
    if not os.path.exists(inputDir):
        os.makedirs(inputDir)


def make_parent_directory(inputFile):
    """
    make the parent directory of the inputFile if it does not exist
    """
    absInputFile = os.path.abspath(inputFile)
    parentDir_absInputFile = os.path.dirname(absInputFile)
    make_directory(parentDir_absInputFile)


def compress_index_bed(inputFile, outputFile, bgzip_cmd, tabix_cmd):

    ####################
    # compress by bgzip
    hOUT = open(outputFile, "w")
    subprocess.call([bgzip_cmd, "-f", outputFile], stdout = hOUT)
    hOUT.close()
    ####################

    ####################
    # index by tabix
    subprocess.call([tabix_cmd, "-p", "bed", outputFile])
    ####################


def sortBedpe(inputFile, outputFile):

    hOUT = open(outputFile, "w")
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", inputFile], stdout = hOUT)
    hOUT.close()

