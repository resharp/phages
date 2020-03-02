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
# 4. sample_measures.txt    (extra generated measures aggregated from results by CalcDiversiMeasures)
class MakeSamplePlots:

    project_dir = ""
    dir_sep = ""
    plot_dir = ""

    sample_meta_df = None
    sample_stats_df = None
    ref_meta_df = None
    sample_measures_df = None
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

        self.read_sample_measures()

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
                                           , names=["run", "ref", "mapped", "nonmapped"])

        # rename NC_024711.1 to crassphage_refseq, for being able to join with sample_measures
        self.sample_stats_df.loc[self.sample_stats_df.ref == "NC_024711.1", "ref"] = "crassphage_refseq"

    def read_ref_metadata(self):

        # NB: ref_genome_ids.txt should be copied from scripts/mgx/ref_seqs
        ref_file_name = self.sample_dir + self.dir_sep + "ref_genome_ids.txt"

        self.ref_meta_df = pd.read_csv(ref_file_name
                                       , sep="\t"
                                       , header=None
                                       , comment="#"
                                       , names=["ref", "genus"]
                                       )

    def read_sample_measures(self):

        # sample measures contains the breadth_1x, ..10x, ..50x and ..100x fractions
        measure_file_name = self.sample_dir + self.dir_sep + "sample_measures.txt"
        self.sample_measures_df = pd.read_csv(measure_file_name,
                                              sep="\t",
                                              )
        nr_1x_05 = len(self.sample_measures_df[self.sample_measures_df.breadth_1x > 0.05])
        nr_1x_20 = len(self.sample_measures_df[self.sample_measures_df.breadth_1x > 0.20])
        nr_1x_80 = len(self.sample_measures_df[self.sample_measures_df.breadth_1x > 0.80])
        nr_1x_95 = len(self.sample_measures_df[self.sample_measures_df.breadth_1x > 0.95])

        logging.info("nr of 1x > 5%: {nr}".format(nr=nr_1x_05))
        logging.info("nr of 1x > 20%: {nr}".format(nr=nr_1x_20))
        logging.info("nr of 1x > 80%: {nr}".format(nr=nr_1x_80))
        logging.info("nr of 1x > 95%: {nr}".format(nr=nr_1x_95))

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

        merge_df = merge_df.merge(self.sample_measures_df
                                  , left_on=[merge_df.run, merge_df.ref]
                                  # nb: sample is reserved word of pandas DataFrame
                                  , right_on=[self.sample_measures_df["sample"], self.sample_measures_df.ref]
                                  , how="left").drop(["key_0", "key_1", "ref_y", "sample"], axis=1)
        merge_df.rename(columns={'ref_x': 'ref'}, inplace=True)

        self.merge_df = merge_df

    def make_bar_plot(self):

        figure_name = "{}{}sample_plots.totals.mapped.genus.pdf".format(self.plot_dir, self.dir_sep)

        sns.barplot(self.genus_sorted_df.genus, self.genus_sorted_df.mapped_sum, log=True
                    , order=self.genus_sorted_df.genus)
        plt.title("total reads mapped for genera")
        plt.savefig(figure_name)
        plt.clf()

    def prepare_data(self):

        self.merge_df = self.merge_df[self.merge_df.mapped != 0]

        # normalize to nr mapped reads per 50 million total reads
        self.merge_df['mapped_norm'] = 50e6 * self.merge_df.mapped / \
                                    (self.merge_df.nonmapped + self.merge_df.mapped)

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
        merge_df["family"] = merge_df.apply(self.family, axis=1)

        merge_df.loc[merge_df.age_cat_short == "B", "age_cat"] = "1. baby"
        merge_df.loc[merge_df.age_cat_short == "4M", "age_cat"] = "2. 4 months"
        merge_df.loc[merge_df.age_cat_short == "12M", "age_cat"] = "3. 12 months"
        merge_df.loc[merge_df.age_cat_short == "M", "age_cat"] = "4. mother"

        self.merge_df = merge_df.sort_values("age_cat")

    @staticmethod
    def family(row):
        return row.sample_name.split("_")[0]

    @staticmethod
    def age_category(row):
        return row.sample_name.split("_")[1]

    @staticmethod
    def shorten_genus(row):
        return row.genus.replace("genus_", "")

    # this is a specific plot for the project that contains metadata about age in the sample name
    def make_abundance_plots(self):

        self.make_category_plot("swarm")
        self.make_category_plot("box")

    def make_category_plot(self, kind):
        sns.catplot(x="genus", y="log10_mapped", kind=kind, data=self.merge_df, hue="age_cat",
                    palette=sns.color_palette("coolwarm", 4), order=self.genus_sorted_df.genus)

        plt.title("log abundance of mapped reads for one ref genome per genus")

        figure_name = "{}{}sample_plots.categories.mapped.genus.{kind}.pdf".format(self.plot_dir, self.dir_sep,
                                                                                   kind=kind)
        plt.savefig(figure_name)
        plt.clf()

    @staticmethod # note: breadth is a pandas series, ge means greater or equal than
    def ge_1_perc(breadth):
        return breadth[breadth.ge(0.01)].count().astype(int)

    @staticmethod
    def ge_5_perc(breadth):
        return breadth[breadth.ge(0.05)].count().astype(int)

    @staticmethod
    def ge_20_perc(breadth):
        return breadth[breadth.ge(0.2)].count().astype(int)

    @staticmethod
    def ge_80_perc(breadth):
        return breadth[breadth.ge(0.8)].count().astype(int)

    @staticmethod
    def ge_95_perc(breadth):
        return breadth[breadth.ge(0.95)].count().astype(int)

    def make_filter_sample_plots(self):
        # now we want to show the number of valid samples after applying filter
        # for each genus we want to show a barplot with nr of samples at different cut-offs
        # so we want to show
        counts_df = self.merge_df.groupby(["ref", "genus"]).agg(
            {
                "breadth_1x": [self.ge_1_perc, self.ge_5_perc, self.ge_20_perc, self.ge_80_perc, self.ge_95_perc],
                "breadth_10x": [self.ge_1_perc, self.ge_5_perc, self.ge_20_perc, self.ge_80_perc, self.ge_95_perc],
                "breadth_50x": [self.ge_1_perc, self.ge_5_perc, self.ge_20_perc, self.ge_80_perc, self.ge_95_perc],
                "breadth_100x": [self.ge_1_perc, self.ge_5_perc, self.ge_20_perc, self.ge_80_perc, self.ge_95_perc]
            }
        ).reset_index()

        counts_df.columns = ["_".join(x) for x in counts_df.columns.ravel()]
        counts_df.rename(columns={'genus_': 'genus', "ref_": "ref"}, inplace=True)

        depths = ["1x", "10x", "50x", "100x"]
        percentages = [1, 5, 20, 80, 95]
        for depth in depths:
            for perc in percentages:
                sns.barplot(x=counts_df.genus,
                            y=counts_df["breadth_{depth}_ge_{perc}_perc".format(depth=depth, perc=perc)])

                figure_name = "{dir}{sep}sample_plots.filter.{depth}.{perc}_perc.pdf".format(
                    dir=self.plot_dir, sep=self.dir_sep, depth=depth, perc=perc
                )
                plt.savefig(figure_name)
                plt.clf()

        debug = True


    def make_scatter_plots(self):
        # total number of mappings between age 4 and 12 months and mother?
        # might also make a mixed scatterplot matrix

        df_families = self.merge_df.pivot_table(
            values="mapped_norm",
            index=["family", "genus"],
            columns="age_cat"
        ).reset_index()

        # extra derived field: hue does not work with strings that only contains numbers
        df_families["genus_"] = "genus_" + df_families["genus"]

        # rename because order does not matter anymore
        df_families.rename(columns={'1. baby': 'baby',
                                    '2. 4 months': '4 months',
                                    '3. 12 months': '12 months',
                                    '4. mother': 'mother'}, inplace=True)

        genera = df_families.genus.unique().tolist()

        for genus in genera:
            self.make_scatter_plot(df_families, genus, "mother", "baby")
            self.make_scatter_plot(df_families, genus, "baby", "4 months")
            self.make_scatter_plot(df_families, genus, "4 months", "12 months")
            self.make_scatter_plot(df_families, genus, "mother", "12 months")

    def make_scatter_plot(self, data, genus, x_value, y_value):

        data = data[data.genus == genus]

        # s_plot = sns.scatterplot(x=x_value, y=y_value, data=data, hue="genus_")
        s_plot = sns.scatterplot(x=x_value, y=y_value, data=data)
        s_plot.set(xscale="log")
        s_plot.set(yscale="log")

        plt.title("genus {genus}: normalized mapped reads/50 million total reads".format(genus=genus))

        y_value = y_value.replace(" ", "_")
        figure_name = "{dir}{sep}sample_plots.scatter.{x_value}.{y_value}.genus_{genus}.png".format(
            dir=self.plot_dir, sep=self.dir_sep,
            x_value=x_value, y_value=y_value, genus=genus
        )
        plt.savefig(figure_name)
        s_plot.clear()
        plt.close(s_plot)

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

        self.make_abundance_plots()

        self.make_filter_sample_plots()

        verbose = False
        # some verbose stuff:
        if verbose:
            self.make_scatter_plots()


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