import argparse
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns;sns.set()
from scipy.stats import entropy
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
    genus_df = None

    genus_palette = {}

    age_cat_palette = {}

    genus_order = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]

    filter = ""

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

        # to do: also read from ref_dir
        # NB: for now ref_genome_ids.txt should be copied from scripts/mgx/ref_seqs
        ref_file_name = self.sample_dir + self.dir_sep + "ref_genome_ids.txt"

        self.ref_meta_df = pd.read_csv(ref_file_name
                                       , sep="\t"
                                       , header=None
                                       , comment="#"
                                       , names=["ref", "genus"]
                                       )
        self.ref_meta_df.genus = self.ref_meta_df.apply(self.shorten_genus, axis=1)

        # to do: cleaner way to exclude this row
        self.ref_meta_df = self.ref_meta_df[self.ref_meta_df.ref != "NC_024711.1"]

        # add color palettes the same way as in MakeGenePlots
        self.build_color_palette_from_ref()

    # use same color coding in MakeGenePlots for consistent color coding of genera and accompanying refs
    def build_color_palette_from_ref(self):
        # If you would like to extend manually on the palette colors and copy color codes see:
        # https://github.com/mwaskom/seaborn/blob/master/seaborn/palettes.py
        colors = sns.color_palette("bright", n_colors=10)
        i = 0
        for index, row in self.ref_meta_df.iterrows():
            self.genus_palette[row.genus] = colors[i]
            i = i + 1

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

        return merge_df

    def filter_on_breadth_threshold(self):

        merge_df = self.merge_df

        breadth_threshold = 0.05
        self.filter = "{perc}%/1x".format(perc=breadth_threshold*100)

        return merge_df[merge_df.breadth_1x > breadth_threshold]

    def make_bar_plot(self):

        figure_name = "{}{}sample_plots.totals.mapped.genus.pdf".format(self.plot_dir, self.dir_sep)

        ax = sns.barplot(self.genus_df.genus, self.genus_df.mapped_sum, log=True
                         , order=self.genus_order
                         , palette=self.genus_palette)
        ax.set(xlabel='crAss-like genus', ylabel='total # of reads mapped')

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

        self.merge_df["log10_mapped"] = np.log10(self.merge_df.mapped_norm)

        self.merge_df["genus"] = self.merge_df.apply(self.shorten_genus, axis=1)

        return self.merge_df

    def sort_genus_according_to_abundance(self):

        df = self.merge_df.groupby(["ref", "genus"]).agg(
            {
                'mapped': ["sum", "std"]
            }).reset_index()
        df.columns = ["_".join(x) for x in df.columns.ravel()]
        df.rename(columns={'genus_': 'genus', "ref_": "ref"}, inplace=True)

        return df

    @staticmethod
    def ge_5_perc(breadth):
        return breadth[breadth.ge(0.05)].count().astype(int)

    @staticmethod
    def ge_50_perc(breadth):
        return breadth[breadth.ge(0.5)].count().astype(int)

    @staticmethod
    def ge_95_perc(breadth):
        return breadth[breadth.ge(0.95)].count().astype(int)

    def make_filter_sample_plots(self):
        # now we want to show the number of valid samples after applying filter
        # for each genus we want to show a barplot with nr of samples at different cut-offs
        # so we want to show
        counts_df = self.merge_df.groupby(["ref", "genus"]).agg(
            {
                "breadth_1x": [self.ge_5_perc, self.ge_50_perc, self.ge_95_perc],
                "breadth_10x": [self.ge_5_perc, self.ge_50_perc, self.ge_95_perc],
                "breadth_50x": [self.ge_5_perc, self.ge_50_perc, self.ge_95_perc]
            }
        ).reset_index()

        counts_df.columns = ["_".join(x) for x in counts_df.columns.ravel()]
        counts_df.rename(columns={'genus_': 'genus', "ref_": "ref"}, inplace=True)

        depths = ["1x", "10x", "50x"]
        percentages = [5, 50, 95]
        for depth in depths:
            for perc in percentages:
                ax = sns.barplot(x=counts_df.genus,
                                 y=counts_df["breadth_{depth}_ge_{perc}_perc".format(depth=depth, perc=perc)],
                                 palette=self.genus_palette,
                                 order=self.genus_order
                                 )

                ax.set(xlabel='crAss-like genus', ylabel='#of samples with breadth > {perc}% for {depth}'.
                       format(depth=depth, perc=perc))
                plt.title("number of samples in which genera are detected")

                figure_name = "{dir}{sep}sample_plots.filter.{depth}.{perc}_perc.pdf".format(
                    dir=self.plot_dir, sep=self.dir_sep, depth=depth, perc=perc
                )
                plt.savefig(figure_name)
                plt.clf()

    # specific ERP005989 enrichment
    # prepare data plot for the project that contains metadata about age in the sample name
    def prepare_specifics_for_project(self):

        merge_df = self.merge_df
        merge_df["age_cat_short"] = merge_df.apply(self.age_category, axis=1)
        merge_df["family"] = merge_df.apply(self.family, axis=1)

        merge_df.loc[merge_df.age_cat_short == "B", "age_cat"] = "baby"
        merge_df.loc[merge_df.age_cat_short == "4M", "age_cat"] = "4 months"
        merge_df.loc[merge_df.age_cat_short == "12M", "age_cat"] = "12 months"
        merge_df.loc[merge_df.age_cat_short == "M", "age_cat"] = "mother"

        merge_df = merge_df.sort_values("age_cat")

        return merge_df

    def write_sample_genera_matrix(self):

        # self.merge_df now contains among others:
        # run (= sample); age_cat (dimension of run/sample)
        # ref, genus
        # value: mapped_norm (normalized)

        df = self.merge_df
        df["ref_genus"] = df["ref"] + " (" + df["genus"] + ")"
        df = df[["run", "ref_genus", "mapped_norm"]]
        df = df.sort_values(["run","ref_genus"])

        df = df.set_index(["run", "ref_genus"])

        # convert from multi-index to matrix
        df = df.unstack()
        df.columns = ["_".join(x) for x in df.columns.ravel()]
        df.columns = [x.replace("mapped_norm_","") for x in df.columns]
        df = df.reset_index()

        df_age = self.merge_df[["run", "age_cat"]]
        df_age = df_age.sort_values("run")
        df_age.reset_index()

        # to do: fix this for adding the age category (somehow it results in double rows (because of multiindex?)
        # df = df_age.merge(df
        #                   , left_on=df_age.run
        #                   , right_on=df.run
        #                   , how="inner").drop(["key_0", "run_y"], axis=1)
        #
        # df.rename(columns={'run_x': 'sample'}, inplace=True)

        file_name = "{}{}sample_genera_matrix.txt".format(
            self.plot_dir, self.dir_sep,
        )
        df.to_csv(path_or_buf=file_name, sep='\t', index=False)

    def build_color_palette_for_ages(self):
        # If you would like to extend manually on the palette colors and copy color codes see:
        # https://github.com/mwaskom/seaborn/blob/master/seaborn/palettes.py
        # check color codes: https://www.color-hex.com/color/00d7ff

        green = "#138D75"
        dark_blue = "#2E86C1"
        purple = "#884EA0"
        orange = "#F39C12"

        bright_4 = [green, dark_blue, purple, orange]

        self.age_cat_palette["baby"] = green
        self.age_cat_palette["4 months"] = dark_blue
        self.age_cat_palette["12 months"] = purple
        self.age_cat_palette["mother"] = orange

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

        # to do: we should filter out the reads that are not above the 5% 1x threshold
        sns.catplot(x="genus",
                    y="log10_mapped",
                    kind=kind,
                    data=self.merge_df,
                    hue="age_cat",
                    hue_order=["baby", "4 months", "12 months", "mother"],
                    palette=self.age_cat_palette, order=self.genus_order
                    )

        plt.title("log abundance of normalized mapped reads for one ref genome per genus")
        plt.ylabel("log10 (norm. mapped reads)")

        figure_name = "{}{}sample_plots.categories.mapped.genus.{kind}.pdf".format(self.plot_dir, self.dir_sep,
                                                                                   kind=kind)
        plt.savefig(figure_name)
        plt.clf()

    @staticmethod
    def macro_entropy(mean_depth):

        macro_entropy = entropy(mean_depth, base=10)

        return macro_entropy

    @staticmethod
    def nr_genera(mean_depth):

        return len(mean_depth)

    def make_diversity_plots(self):

        # to do: use genus and mean_depth to calculate entropy for each sample and categorize for age_cat
        df_entropies = self.merge_df.groupby(["run", "age_cat", "age_cat_short" ]).agg(
            {
                'mean_depth': ["mean", self.nr_genera, self.macro_entropy]
            }
        ).reset_index()

        df_entropies.columns = ["_".join(x) for x in df_entropies.columns.ravel()]
        df_entropies.rename(columns={'run_': 'run',
                                     'age_cat_': 'age_cat',
                                     "age_cat_short_": "age_cat_short",
                                     "mean_depth_nr_genera": "nr_genera"
                                     }, inplace=True)

        self.make_cat_plot_for_age_categories(df_entropies, "mean_depth_macro_entropy",
                                              "macro diversity (Shannon entropy)")
        self.make_cat_plot_for_age_categories(df_entropies, "nr_genera",
                                              "macro diversity (nr of genera)")

        measure="nr_genera"
        title="macro diversity: nr of genera"
        self.regression_plot(df_entropies, measure, title)

        measure="mean_depth_macro_entropy"
        title="macro diversity: entropy based on genus abundances"
        self.regression_plot(df_entropies, measure, title)

    def make_cat_plot_for_age_categories(self, data, measure, title):

        kind = "swarm"
        sns.catplot(x="age_cat", y=measure, kind=kind, data=data,
                    palette=self.age_cat_palette
                    , order=["baby", "4 months", "12 months", "mother"]
                    )

        plt.title("{title} ({filter})".format(
            title=title, filter=self.filter
        ))

        figure_name = "{dir}{sep}catplot.age_cat.{measure}.{filter}.svg".format(
            dir=self.plot_dir, sep=self.dir_sep, measure=measure, filter=self.filter.replace("/", "_"),
        )
        plt.savefig(figure_name)
        plt.clf()

    def regression_plot(self, data, measure, title):

        data = data[data.age_cat_short != "B"]

        g0 = sns.lmplot(x="mean_depth_mean", y=measure,
                        palette=self.age_cat_palette,
                        hue="age_cat",
                        hue_order=["4 months", "12 months", "mother"],
                        data=data,
                        height=5,
                        aspect=1.5)

        if measure == "nr_genera":
            plt.ylim(0, 7)
        if measure == "mean_depth_macro_entropy":
            plt.ylim(-0.2, 0.8)
            plt.ylabel("macro entropy")

        title = "{title} ({filter})".format(title=title, filter=self.filter)
        plt.xlabel("mean sample depth")
        plt.title(title)

        figure_name = "{}{}macro_{measure}_against_depth_{filter}.svg".\
            format(self.plot_dir, self.dir_sep, measure=measure, filter=self.filter.replace("/", "_"))

        plt.savefig(figure_name)
        plt.clf()

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

        genera = df_families.genus.unique().tolist()

        for genus in genera:
            self.make_scatter_plot(df_families, genus, "mother", "baby")
            self.make_scatter_plot(df_families, genus, "baby", "4 months")
            self.make_scatter_plot(df_families, genus, "4 months", "12 months")
            self.make_scatter_plot(df_families, genus, "mother", "12 months")

    def make_scatter_plot(self, data, genus, x_value, y_value):

        data = data[data.genus == genus]

        # s_plot = sns.scatterplot(x=x_value, y=y_value, data=data, hue="genus_")
        s_plot = sns.scatterplot(x=x_value, y=y_value, data=data, palette=self.age_cat_palette)
        s_plot.set(xscale="log")
        s_plot.set(yscale="log")

        plt.title("genus {genus}: normalized mapped reads/50 million total reads".format(genus=genus))

        y_value = y_value.replace(" ", "_")
        figure_name = "{dir}{sep}sample_plots.scatter.{x_value}.{y_value}.genus_{genus}.png".format(
            dir=self.plot_dir, sep=self.dir_sep,
            x_value=x_value, y_value=y_value, genus=genus
        )
        plt.savefig(figure_name)
        plt.clf()

    def create_plot_dir(self):

        self.plot_dir = self.sample_dir + self.dir_sep + "SamplePlots"
        os.makedirs(self.plot_dir, exist_ok=True)

    def do_analysis(self):

        self.read_files()

        self.create_plot_dir()

        self.merge_df = self.merge_files()

        self.merge_df = self.prepare_data()

        self.genus_df = self.sort_genus_according_to_abundance()

        self.build_color_palette_for_ages()

        self.make_filter_sample_plots()

        self.merge_df = self.filter_on_breadth_threshold()

        self.make_bar_plot()

        self.merge_df = self.prepare_specifics_for_project()

        self.write_sample_genera_matrix()

        self.make_abundance_plots()

        self.make_diversity_plots()

        verbose = False

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

# to do for testing, do not use in production
project_dir= r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"
do_analysis(["-d", project_dir])
