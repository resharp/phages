import argparse
import logging
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns;sns.set()
from scipy.stats import entropy

# We will calculate per protein position (by pile up/stacking of all reads of all samples on top of each other):
# the overall coverage,
# the total number of possible codon types for that position (may be 1)
# and the entropy at that position (over all samples), [not including the ref genome codon type!]
# (we can later aggregate this per gene, to have less data points).
# Goal is to have a trend line for the dependence of the entropy or SNP density (micro diversity)
#
# results will be put in CodonMeasures directory (sub dir of sample_dir)
# plots will be made by different class

# Also, we want to split the samples based on age_category, and therefore we need the sample metadata.
# First calculate the overall entropy, and later make a stack for the age categories.

class CalcCodonMeasures:

    sample_dir = ""
    ref = ""

    aa_df = None
    sample_meta_df = None
    sample_measures_df = None

    def __init__(self, sample_dir, ref):

        self.sample_dir = sample_dir
        self.ref = ref

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'CalcCodonMeasures.log', filemode='w',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level=logging.DEBUG)

    def read_and_concat_measures(self, file_name_ext, ref=None):

        subfolders = [f.path for f in os.scandir(self.sample_dir) if f.is_dir() and ref in f.name]

        samples = []

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            sample = sample.split("_")[0]
            sample_name = subfolder + self.dir_sep + sample + file_name_ext

            # to do: check if the sample is above breadth threshold
            # based on sample measures

            if os.path.isfile(sample_name) and os.path.getsize(sample_name) > 0:
                logging.info("processing {}".format(sample_name))

                measure_df = pd.read_csv(sample_name
                                         ,  sep='\t'
                                         )
                samples.append(measure_df)

        return pd.concat(samples)

    def read_files(self):

        logging.debug("start reading tables")
        # to do: read sample measures (for filter criteria, e.g. 5%/10x)

        # to do: read sample metadata (for age categories)

        # to do: read AA measures

        logging.debug("end reading tables")

        pass

    def calc_measures(self):
        pass

    def create_output_dir(self):
        pass

    def write_measures(self):
        pass

    def run_calc(self):

        self.read_files()

        self.calc_measures()

        self.create_output_dir()

        self.write_measures()


def run_calc(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    parser.add_argument("-r", "--ref", dest="ref",
                        help="reference genome id", metavar="[ref}", required=True)

    args = parser.parse_args(args_in)

    print("Start running CalcCodonMeasures")
    print("sample_dir: " + args.sample_dir)
    print("reference genome: " + args.ref)

    calc = CalcCodonMeasures(args.sample_dir, args.ref)

    calc.run_calc()

# if __name__ == "__main__":
#     run_calc(sys.argv[1:])


# TODO for testing, do not use in production
sample_dir = r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"

# or run one sample, or a list of
ref = "crassphage_refseq"
run_calc(["-d", sample_dir, "-r", ref])