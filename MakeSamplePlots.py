import argparse
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns;sns.set()
import sys

# MakeSamplePlots
# goal: show sample statistics. How many reads are mapped to what reference genomes?
# how does it relate to some metadata?
# this should be generic for the metadata file format
# however, we might want to plug in some specific stuff for certain projects
# input:
# 1. sample_stats.txt       (statistics of samples vs ref genomes)
# 2. metadata_ERP005989.txt (metadata project)
# 3. ref_genome_ids.txt     (metadata ref genomes)
# 4. sample_ref_measures.txt(extra generated measures aggregated from results by CalcDiversiMeasures)


class MakeSamplePlots:

    project_dir = ""
    dir_sep = ""
    plot_dir = ""

    sample_meta_df = None
    sample_stats_df = None
    ref_meta_df = None
    merge_df = None
    genus_sorted_df = None

    def __init__(self, sample_dir):

        self.sample_dir = sample_dir
        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'MakeSamplePlots.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.INFO)

    def read_files(self):

        logging.info("start reading tables")

        self.read_sample_metadata()

        self.read_sample_stats()

        self.read_ref_metadata()

        logging.info("finished reading tables")

    def read_sample_metadata(self):

        metadata_filename = self.sample_dir+ self.dir_sep + "metadata_ERP005989.txt"

        self.sample_meta_df = pd.read_csv(metadata_filename
                                          , sep='\t'
                                          )

        self.sample_meta_df = self.sample_meta_df[['analysis.run.accession', 'sample.sample_name', 'sample.sample_desc']]

        self.sample_meta_df.rename(columns={'analysis.run.accession': 'run',
                                            'sample.sample_name': 'sample_name',
                                            'sample.sample_desc': 'sample_desc'}, inplace=True)

    def read_sample_stats(self):

        stats_file_name = self.sample_dir + self.dir_sep + "sample_stats.txt"

        self.sample_stats_df = pd.read_csv(stats_file_name
                                           , sep='\t'
                                           , header=None
                                           , names=["run", "ref", "mapped", "nonmapped"]
                                           )

    def read_ref_metadata(self):

        # NB: ref_genome_ids.txt should be copied from scripts/mgx/ref_seqs
        ref_file_name = self.sample_dir + self.dir_sep + "ref_genome_ids.txt"

        self.ref_meta_df = pd.read_csv(ref_file_name
                                       , sep="\t"
                                       , header=None
                                       , comment="#"
                                       , names=["ref", "genus"]
                                       )

    def merge_files(self):

        # now we join all the files together, and this is still generic i.e. not specific for a certain project
        # and after doing that we do the specific ERP005989 enrichment
        merge_df = self.sample_stats_df.merge(self.sample_meta_df
                                              ,   left_on=self.sample_stats_df.run
                                              ,   right_on=self.sample_meta_df.run
                                              ,   how="inner").drop(["key_0", "run_y"], axis=1)
        merge_df.rename(columns={'run_x': 'run'}, inplace=True)

        assert(len(self.sample_stats_df) == len(merge_df))

        merge_df = merge_df.merge(self.ref_meta_df
                                  , left_on=merge_df.ref
                                  , right_on=self.ref_meta_df.ref
                                  , how="inner").drop(["key_0", "ref_y"], axis=1)
        merge_df.rename(columns={'ref_x': 'ref'}, inplace=True)
        assert(len(self.sample_stats_df) == len(merge_df))

        nr_of_mappings = len(merge_df)
        logging.info("Number of mappings analyzed: {nr_of_mappings}".format(nr_of_mappings=nr_of_mappings))
        logging.info("total number of mapped reads: {mapped}".format(mapped=merge_df.mapped.sum()))

        self.merge_df = merge_df

    def make_bar_plot(self):

        figure_name = "{}{}sample_plots.totals.mapped.genus.pdf".format(self.plot_dir, self.dir_sep)

        sns.barplot(self.genus_sorted_df.genus, self.genus_sorted_df.mapped_sum, log=True
                    , order=self.genus_sorted_df.genus)
        plt.title("total reads mapped for genera")
        plt.savefig(figure_name)
        plt.clf()

    def prepare_data(self):

        # merge_df = self.merge_df
        self.merge_df = self.merge_df[self.merge_df.mapped != 0]

        nr_mappings = len(self.merge_df)
        logging.info("nr of non-zero mappings among samples: {nr_mappings}".format(nr_mappings=nr_mappings))

        self.merge_df["log10_mapped"] = np.log10(self.merge_df.mapped)

        self.merge_df["genus"] = self.merge_df.apply(self.shorten_genus, axis=1)

        self.sort_genus_according_to_abundance()

    def sort_genus_according_to_abundance(self):

        df = self.merge_df.groupby(["ref", "genus"]).agg(
            {
                'mapped': ["sum", "std"]
            }).reset_index()
        df.columns = ["_".join(x) for x in df.columns.ravel()]
        df.rename(columns={'genus_': 'genus', "ref_": "ref"}, inplace=True)

        self.genus_sorted_df = df.sort_values('mapped_sum', ascending=False)

        debug = "True"

    # specific ERP005989 enrichment
    # prepare data plot for the project that contains metadata about age in the sample name
    def prepare_specifics_for_project(self):

        merge_df = self.merge_df
        merge_df["age_cat_short"] = merge_df.apply(self.age_category, axis=1)

        merge_df.loc[merge_df.age_cat_short == "B", "age_cat"] = "1. baby"
        merge_df.loc[merge_df.age_cat_short == "4M", "age_cat"] = "2. 4 months"
        merge_df.loc[merge_df.age_cat_short == "12M", "age_cat"] = "3. 12 months"
        merge_df.loc[merge_df.age_cat_short == "M", "age_cat"] = "4. mother"

        self.merge_df = merge_df.sort_values("age_cat")

    @staticmethod
    def age_category(row):
        return row.sample_name.split("_")[1]

    @staticmethod
    def shorten_genus(row):
        return row.genus.replace("genus_", "")

    # this is a specific plot for the project that contains metadata about age in the sample name
    def make_category_plot(self):

        sns.catplot(x="genus", y="log10_mapped", kind="swarm", data=self.merge_df, hue="age_cat",
                    palette=sns.color_palette("coolwarm", 4), order=self.genus_sorted_df.genus)

        plt.title("log abundance of mapped reads for one ref genome per genus")

        figure_name = "{}{}sample_plots.categories.mapped.genus.pdf".format(self.plot_dir, self.dir_sep)
        plt.savefig(figure_name)
        plt.clf()

        sns.catplot(x="genus", y="log10_mapped", kind="box", data=self.merge_df, hue="age_cat",
                    palette=sns.color_palette("coolwarm", 4), order=self.genus_sorted_df.genus)

        plt.title("log abundance of mapped reads for one ref genome per genus")

        figure_name = "{}{}sample_plots.box.categories.mapped.genus.pdf".format(self.plot_dir, self.dir_sep)
        plt.savefig(figure_name)
        plt.clf()

    def create_plot_dir(self):

        self.plot_dir = self.sample_dir + self.dir_sep + "SamplePlots"
        os.makedirs(self.plot_dir, exist_ok=True)

    def do_analysis(self):

        self.read_files()

        self.create_plot_dir()

        self.merge_files()

        self.prepare_data()

        self.make_bar_plot()

        self.prepare_specifics_for_project()

        self.make_category_plot()


def do_analysis(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--project_dir", dest="project_dir",
                        help="project directory with sample vs ref runs in subfolders", metavar="[project_dir]", required=True)

    args = parser.parse_args(args_in)

    print("Start running MakeSamplePlots")
    print("project_dir: " + args.project_dir)

    make = MakeSamplePlots(args.project_dir)

    make.do_analysis()


# if __name__ == "__main__":
#     do_analysis(sys.argv[1:])

#TODO for testing, do not use in production
project_dir= r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"
do_analysis(["-d", project_dir])

# D:\17 Dutihl Lab\_tools\_pipeline\ERP005989