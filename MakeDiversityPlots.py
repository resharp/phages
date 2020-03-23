import argparse
import logging
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns;sns.set()
from scipy.stats import linregress


# MakeDiversityPlots
# 1. plot entropy against average coverage for genes and for all AA positions over all samples
#   for different age categories
# 2. determine significance of trend lines (and significance of distributions?)
class MakeDiversityPlots:

    sample_dir = ""
    ref = ""
    plot_dir = ""

    aa_df = None
    protein_df = None

    def __init__(self, sample_dir, ref):

        self.sample_dir = sample_dir
        self.ref = ref

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        self.plot_dir = self.sample_dir + self.dir_sep + "CodonMeasures"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'MakeDiversityPlots.log', filemode='w',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level=logging.DEBUG)

    def read_and_concat_measures(self, file_name_prefix, ref=None):

        files = [f.path for f in os.scandir(self.plot_dir)
                 if ref in f.name
                 and file_name_prefix in f.name]
        age_sets = []
        age_cats = []

        for file in files:
            base_name = os.path.basename(file)
            age_cat = base_name.split("_")[-1].replace(".txt", "")
            age_cats.append(age_cat)

            logging.info("processing {}".format(file))

            measure_df = pd.read_csv(file
                                     ,  sep='\t'
                                     )
            measure_df["age_cat"] = age_cat
            age_sets.append(measure_df)

        return pd.concat(age_sets, keys=age_cats)

    def read_files(self):

        logging.debug("start reading tables")

        self.aa_df = self.read_and_concat_measures("codon_entropy", self.ref)

        logging.debug("end reading tables")

    def create_plot_dir(self):

        os.makedirs(self.plot_dir, exist_ok=True)

    def aggregate_on_protein_level(self):

        self.protein_df = self.aa_df.groupby(["age_cat", "protein"]).agg(
            {'coverage': 'mean',
             'entropy': 'mean'}).reset_index()

    def make_plots(self):

        title = "crAssphage: mean entropy of all position against mean coverage for that position."
        self.plot_entropy_for_age_categories(self.aa_df, "codon", title)
        title = "crAssphage: mean entropy of all genes against mean coverage for that gene."
        self.plot_entropy_for_age_categories(self.protein_df, "protein", title)

    def plot_entropy_for_age_categories(self, data, level, title):

        data = data[data.age_cat != "all"]

        plt.figure(figsize=(14, 10))
        g0 = sns.lmplot(x="coverage", y="entropy",
                        hue="age_cat",
                        data=data)
        plt.title(title)

        filename="{}{}entropy_against_coverage_for_{level}_{ref}.svg".format(
            self.plot_dir, self.dir_sep, level=level, ref=self.ref)

        plt.savefig(filename)
        plt.clf()

    def write_measures(self):
        # https://stackoverflow.com/questions/33490833/display-regression-equation-in-seaborn-regplot

        # to do: write measures for the slope of the lines for different age categories
        # and determine if they are significantly different
        self.write_measures_for_level(self.protein_df, "protein")
        self.write_measures_for_level(self.aa_df, "codon")

    def write_measures_for_level(self, data, level):

        df_lr = pd.DataFrame(columns=('ref', 'slope', 'intercept', 'r_value', 'p_value', 'std_err'))
        df_lr.index.name = 'age_cat'

        age_cats = data.age_cat.unique()

        for age_cat in age_cats:

            data_age = data[data.age_cat == age_cat]

            x_data = data_age.coverage
            y_data = data_age.entropy

            slope, intercept, r_value, p_value, std_err = linregress(x=x_data, y=y_data)

            df_lr.loc[age_cat] = [
                self.ref, slope, intercept, r_value, p_value, std_err
            ]

        filename = "{}{}linear_regression_for_{level}_{ref}.txt".format(
            self.plot_dir, self.dir_sep, level=level, ref=ref)
        df_lr.to_csv(path_or_buf=filename, sep='\t')

    def do_analysis(self):

        self.read_files()

        self.aggregate_on_protein_level()

        self.create_plot_dir()

        self.make_plots()

        self.write_measures()


def do_analysis(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    parser.add_argument("-r", "--ref", dest="ref",
                        help="reference genome id", metavar="[ref}", required=True)

    args = parser.parse_args(args_in)

    print("Start running MakeDiversityPlots")
    print("sample_dir: " + args.sample_dir)
    print("reference genome: " + args.ref)

    make = MakeDiversityPlots(args.sample_dir, args.ref)

    make.do_analysis()

# if __name__ == "__main__":
#     run_calc(sys.argv[1:])


# TODO for testing, do not use in production
sample_dir = r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"

# or run one sample, or a list of
ref = "crassphage_refseq"
do_analysis(["-d", sample_dir, "-r", ref])

# refs = ["crassphage_refseq", "sib1_ms_5", "err975045_s_1", "inf125_s_2", "srr4295175_ms_5"]
# for ref in refs:
#     do_analysis(["-d", sample_dir, "-r", ref])