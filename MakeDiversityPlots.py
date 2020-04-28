import argparse
import logging
import numpy as np
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns;sns.set()
from scipy.stats import linregress


# MakeDiversityPlots
# 1. plot entropy against average coverage for genes and for all AA positions over all samples
#   for different age categories
# 2. plot SNP density against average coverage for genes and for all AA positions over all samples
#
# 3. determine significance of trend lines (and significance of distributions?)
#
#
# to do: add gene filter 50% or 95% 10x
# to do: add annotation to self.protein_df (so after aggregation per Protein and age category)
#        so we can compare if measures differ more per age than per region (functional category)
# to do: pN/pS calculation for genes per age category per region
#        requires codon bias table
# to do: violin plot with x=region and hue=age_cat for entropy and/or SNP density
# to do: violin plot for pN/pS


class MakeDiversityPlots:

    sample_dir = ""
    ref_dir = ""
    ref = ""
    plot_dir = ""

    aa_df = None
    protein_df = None
    gene_anno_df = None

    def __init__(self, sample_dir, ref_dir, ref):

        self.sample_dir = sample_dir
        self.ref_dir = ref_dir
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

        if ref == "all":
            files = [f.path for f in os.scandir(self.plot_dir)
                     if file_name_prefix in f.name
                     and "xls" not in f.name]
        else:
            files = [f.path for f in os.scandir(self.plot_dir)
                     if ref in f.name
                     and "xls" not in f.name
                     and file_name_prefix in f.name]

        age_sets = []

        for file in files:
            base_name = os.path.basename(file)
            age_cat = base_name.split("_")[-1].replace(".txt", "")
            ref = base_name.replace(file_name_prefix + "_", "").replace("_" + age_cat, "").replace(".txt", "")

            logging.info("processing {}".format(file))

            try:
                measure_df = pd.read_csv(file
                                         , sep='\t'
                                         )
            except:
                logging.error("Could not read file {file}".format(file=file))
                raise

            measure_df["age_cat"] = age_cat
            measure_df["ref"] = ref

            age_sets.append(measure_df)

        return pd.concat(age_sets)

    def read_gene_annotation(self, ref):

        gene_anno_file_name = self.ref_dir + self.dir_sep + "{ref}_gene_list.txt".format(ref=ref)

        anno_df = pd.read_csv(gene_anno_file_name
                              ,   sep='\t'
                              ,   header=None
                              ,   usecols=[0, 1, 2, 3, 10]
                              ,   names=["Protein", "gene_fam", "gene_fam_annot", "region", "Annotation"]
                              )
        anno_df["ref"] = ref

        return anno_df

    def read_files(self):

        logging.debug("start reading tables")

        self.aa_df = self.read_and_concat_measures("codon_entropy", self.ref)

        self.gene_anno_df = self.read_gene_annotation(self.ref)

        logging.debug("end reading tables")

    def merge_files(self):

        data = self.protein_df

        # add annotation, e.g. region
        merge_df = data.merge(self.gene_anno_df
                              , left_on=data.protein
                              , right_on=self.gene_anno_df.Protein
                              , how='left').drop(["key_0", "Protein", "ref_y"], axis=1)
        merge_df.rename(columns={"ref_x": "ref"}, inplace=True)

        self.protein_df = merge_df

    def create_plot_dir(self):

        os.makedirs(self.plot_dir, exist_ok=True)

    def aggregate_on_protein_level(self):

        # calculate derived CntSnp = CntSyn + CntNonSyn
        self.aa_df["CntSnp"] = self.aa_df["CntSyn"].replace(np.nan, 0) + self.aa_df["CntNonSyn"].replace(np.nan, 0)

        count_snp_median = self.aa_df.CntSnp[self.aa_df["CntSnp"] != 0].median()
        snp_pseudo_count = np.sqrt(count_snp_median) / 2

        self.protein_df = self.aa_df.groupby(["age_cat", "protein", "ref"]).agg(
            {'coverage': 'mean',
             'entropy': 'mean',
             'snp': 'sum',
             'position': 'count',
             'syn': 'sum',
             'non_syn': 'sum',
             'CntNonSyn': 'sum',
             'CntSyn': 'sum',
             }
        ).reset_index()

        self.protein_df["snp_density"] = self.protein_df.snp / self.protein_df.position

        self.protein_df["syn_ratio"] = self.protein_df["syn"] / self.protein_df["non_syn"]

        self.protein_df["pN_pS"] = self.protein_df["syn_ratio"] * \
             (self.protein_df["CntNonSyn"] + snp_pseudo_count)/(self.protein_df["CntSyn"] + snp_pseudo_count)

        self.protein_df["log10_pN_pS"] = np.where(self.protein_df["pN_pS"] > 0, np.log10(self.protein_df["pN_pS"]), 0)

    def make_plots(self):

        # you may want to look at the plot with data points for all AA positions but it is very messy
        # title = "{ref}: mean entropy of all position against mean coverage for that position.".format(ref=self.ref)
        # self.plot_entropy_for_age_categories(self.aa_df, "snp_density", "codon", title)

        title = "{ref}: mean SNP density of all genes against mean coverage for that gene.".format(ref=self.ref)
        self.plot_measure_for_age_categories(self.protein_df, "snp_density", "protein", title)

        title = "{ref}: mean entropy of all genes against mean coverage for that gene.".format(ref=self.ref)
        self.plot_measure_for_age_categories(self.protein_df, "entropy", "protein", title)

    def make_violin_plots(self):
        self.make_violin_plot("log10_pN_pS", "log10(pN/pS)")
        self.make_violin_plot("entropy", "entropy")
        self.make_violin_plot("snp_density", "SNP density")

    def make_violin_plot(self, measure, measure_label):
        # you may want to download the new results from the server
        data = self.protein_df
        data.loc[data.region == "assembly", "region"] = "assembly.rest"

        hue_order = ["B", "4M", "12M", "M", "all"]
        if measure != "log10_pN_pS":
            hue_order.remove("B")   # for SNP density and entropy we remove the B (baby): not enough data

        plt.figure(figsize=(12, 6))
        ax = sns.violinplot(x="region", y=measure,
                            hue="age_cat",
                            data=data, palette="muted",
                            inner="stick",
                            order=["replication", "transcription", "assembly.capsid", "assembly.tail", "assembly.rest"],
                            hue_order=hue_order
                            )
        title = self.ref + ": every stick is info from mapped reads from all samples in age category for one gene"

        ax.set(xlabel='functional region', ylabel=measure_label)

        plt.title(title)

        filename="{}{}{measure}_violin_plots_{ref}.svg".format(
            self.plot_dir, self.dir_sep, measure=measure, ref=self.ref)

        # plt.show()
        plt.savefig(filename)
        plt.clf()

    def plot_measure_for_age_categories(self, data, measure, level, title):

        data = data[data.age_cat != "all"]

        # we exclude the age category "B" because they have a very low coverage
        # and the trend line is messing up the figure
        data = data[data.age_cat != "B"]
        # data = data[data.coverage < 100]
        # data = data[data.ref != "crassphage_refseq"]

        # you can filter out data points below a certain coverage:
        # data = data[data.coverage > 250]

        g0 = sns.lmplot(x="coverage", y=measure,
                        hue="age_cat",
                        data=data,
                        height=5, aspect=1.5)
        if self.ref == "crassphage_refseq" or self.ref == "all":
            # we want these pictures to be in the same format for comparing
            plt.xlim(-1000, 13000)
        if measure == "entropy":
            plt.ylim(-0.01, 0.15)
        if measure == "snp_density":
            plt.ylim(-0.02, 0.80)
        # you may want to look at the residual plots before deciding on any possible trend
        # sns.residplot(x="coverage", y="entropy", data=data[data.age_cat == "M"])

        plt.title(title)

        filename="{}{}{measure}_against_coverage_for_{level}_{ref}.svg".format(
            self.plot_dir, self.dir_sep, measure=measure, level=level, ref=self.ref)

        # plt.show()
        plt.savefig(filename)
        plt.clf()

    def write_measures(self):
        # https://stackoverflow.com/questions/33490833/display-regression-equation-in-seaborn-regplot

        # write measures for the slope of the lines for different age categories
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
            self.plot_dir, self.dir_sep, level=level, ref=self.ref)
        df_lr.to_csv(path_or_buf=filename, sep='\t')

    def do_analysis(self):

        self.read_files()

        self.aggregate_on_protein_level()

        # we merge after aggregation because we only need annotation on Protein level
        self.merge_files()

        self.create_plot_dir()

        self.make_plots()

        self.make_violin_plots()

        self.write_measures()


def do_analysis(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    parser.add_argument("-rd", "--ref_dir", dest="ref_dir",
                        help="directory reference genomes", metavar="[ref_dir]", required=True)

    parser.add_argument("-r", "--ref", dest="ref",
                        help="reference genome id", metavar="[ref]", required=False)

    parser.add_argument("-a", "--all", action="store_true", default=False,
                        help="run analysis on all ref genomes at once", required=False)

    args = parser.parse_args(args_in)

    print("Start running MakeDiversityPlots")
    print("sample_dir: " + args.sample_dir)
    if args.ref:
        print("reference genome: " + args.ref)
    if args.all:
        print("run analysis on all reference genomes")
        args.ref = "all"

    make = MakeDiversityPlots(args.sample_dir, args.ref_dir, args.ref)

    make.do_analysis()



# if __name__ == "__main__":
#     run_calc(sys.argv[1:])


# TODO for testing, do not use in production
sample_dir = r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"
ref_dir = r"D:\17 Dutihl Lab\source\phages_scripts\mgx\ref_seqs"

# # or run one sample, or a list of
ref = "crassphage_refseq"
# ref = "hvcf_a6_ms_4"
# ref = "sib1_ms_5"

do_analysis(["-d", sample_dir, "-rd", ref_dir, "-r", ref])

# refs = ["crassphage_refseq", "sib1_ms_5", "err975045_s_1", "inf125_s_2", "srr4295175_ms_5",
#         "hvcf_a6_ms_4", "fferm_ms_11", "err844030_ms_1", "eld241-t0_s_1", "cs_ms_21"]
# for ref in refs:
#     do_analysis(["-d", sample_dir, "-r", ref])

# do_analysis(["-d", sample_dir, "-a"])