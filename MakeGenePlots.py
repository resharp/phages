import argparse
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns;sns.set()
import sys

# MakeGenePlots
# create multiple heatmaps
#   1st with the log pN/pS values
#   2nd with the missing genes (just 0/1 for clarity)
#   3rd distribution plots [integration over structural genes or other categories]
#
#   all sample_gene measures are in small separate files
#   we aggregate them to gene level (aggregation over all samples)
#   therefore we load all files and merge them into one dataframe
#
# masking missing values in heatmaps:
# https://github.com/mwaskom/seaborn/issues/375


class MakeGenePlots:

    sample_dir = ""
    ref_dir = ""
    dir_sep = ""
    plot_dir = ""

    gene_sample_df = None
    gene_df = None
    sample_df = None
    gene_anno_df = None

    bin_sample_df = None
    bin_df = None

    sample_breadth_df = None

    threshold_depth = 0
    threshold_breadth = 0

    def __init__(self, sample_dir, ref_dir, threshold_depth, threshold_breadth):

        self.sample_dir = sample_dir
        self.ref_dir = ref_dir
        self.threshold_depth = threshold_depth
        self.threshold_breadth = threshold_breadth

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'MakeGenePlots.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.INFO)

    def read_files(self, ref):

        logging.info("start reading tables")

        self.gene_sample_df = self.read_and_concat_sample_measures(ref, "_gene_measures.txt")
        self.read_gene_annotation(ref)

        self.bin_sample_df = self.read_and_concat_sample_measures(ref, "_bin_measures.txt")

        self.sample_breadth_df = self.read_and_concat_sample_measures(ref, "_sample_measures.txt")

        logging.info("finished reading tables")

    def read_and_concat_sample_measures(self, ref, file_name_ext):

        subfolders = [f.path for f in os.scandir(self.sample_dir) if f.is_dir() and ref in f.name]
        samples = []

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            sample = sample.split("_")[0]
            sample_name = subfolder + self.dir_sep + sample + file_name_ext

            if os.path.isfile(sample_name) and os.path.getsize(sample_name) > 0:
                logging.info("processing {}".format(sample_name))

                gene_sample_df = pd.read_csv(sample_name
                                         ,  sep='\t'
                                         )
                samples.append(gene_sample_df)

        return pd.concat(samples)

    def read_gene_annotation(self, ref):

        gene_anno_file_name = self.ref_dir + self.dir_sep + "{ref}_gene_list.txt".format(ref=ref)

        self.gene_anno_df = pd.read_csv(gene_anno_file_name
                                        ,   sep='\t'
                                        ,   header=None
                                        ,   usecols=[0,1]
                                        ,   names=["Protein", "Annotation"]
                                        )
        pass

    def prepare_data_for_plots(self):

        # here we can prepare the data further for displaying purposes
        # it could also have been done in CalcDiversiMeasures

        # we want to see the missing genes, and missing genes are if there are no mappings or hardly any
        self.gene_sample_df['missing_gene'] = 0
        self.gene_sample_df.loc[self.gene_sample_df.AAcoverage_perc < 0.02 , 'missing_gene'] = 1
        self.gene_sample_df.loc[np.isnan(self.gene_sample_df.AAcoverage_perc), 'missing_gene'] = 1

        self.gene_sample_df['double_coverage'] = 0
        self.gene_sample_df.loc[self.gene_sample_df.AAcoverage_perc > 2 , 'double_coverage'] = 1

        # data_debug = self.gene_df[['sample','Protein', 'AAcoverage_perc', 'missing_gene', 'double_coverage']]

        # make binary plot for positive selection or conservation
        self.gene_sample_df['positive_selection'] = 0
        self.gene_sample_df.loc[self.gene_sample_df["log10_pN/pS"] > -0.3, 'positive_selection'] = 1
        self.gene_sample_df.loc[self.gene_sample_df["log10_pN/pS"] < -0.7, 'positive_selection'] = -1

    def create_plot_dir(self, ref):

        self.plot_dir = self.sample_dir + self.dir_sep + ref + "_" + "GenePlots"
        os.makedirs(self.plot_dir, exist_ok=True)

    def create_histograms(self):

        # create a histogram of the occurence of genes across samples
        data = self.gene_sample_filtered_on_quality()

        # we only need the counts per gene

        counts_df = data.groupby("Protein").agg(
            {
                'log10_pN/pS': ["count"]
            }).reset_index()
        counts_df.columns = ["_".join(x) for x in counts_df.columns.ravel()]
        counts_df.rename(columns={'Protein_': 'Protein'}, inplace=True)

        min = counts_df['log10_pN/pS_count'].min()
        max = counts_df['log10_pN/pS_count'].max()

        sns.distplot(counts_df[['log10_pN/pS_count']], bins=(max-min+1), kde=False)

        plt.title("Distribution of sample presence for genes ({breadth}/{depth}x)".
                  format(depth=self.threshold_depth, breadth=self.threshold_breadth))
        figure_name = "{}{}gene_plots.gene_sample_counts.{breadth}.{depth}x.pdf".\
            format(self.plot_dir, self.dir_sep, depth=self.threshold_depth, breadth=self.threshold_breadth)

        plt.savefig(figure_name)

    def create_heatmaps(self):

        # filter gene sample on quality
        data = self.gene_sample_filtered_on_quality()

        # also filter out complete samples before showing pN/pS
        # we do not care for samples with a 1x horizontal coverage below 5%
        filtered_data = self.gene_sample_filtered_on_sample_coverage(data)

        # make a heatmap of the log10_pN/pS based on multiple samples
        self.create_heatmap(filtered_data, "log10_pN/pS", "Log 10 of pN/pS (red=diverging, blue=conserved)")
        self.create_heatmap(filtered_data, "positive_selection", "log10_pN/pS either > -0.3 or < - 0.7")

        self.create_heatmap(filtered_data, "SndAAcnt_perc_polymorphism_mean", "Within sample AA variation in genes")

        self.create_heatmap(filtered_data, "entropy_mean", "Mean codon entropy (base 10)")

        #make a heatmap of quality measure AAcoverage_cv
        self.create_heatmap(data, "AAcoverage_cv", "Internal coefficient of variation per gene")
        self.create_heatmap(data, "AAcoverage_perc", "Coverage percentage (compared to whole genome)")

        #use unfiltered data for quality measures
        data = self.gene_sample_df

        self.create_heatmap(data, "breadth_1x", "1x horizontal coverage percentages")
        self.create_heatmap(data, "breadth_10x", "10x horizontal coverage percentages")
        self.create_heatmap(data, "breadth_50x", "50x horizontal coverage percentages")

        self.create_heatmap(data, "missing_gene", "Missing genes (based on coverage < 2% rest genome)")

        self.create_heatmap(data, "double_coverage", "Genes that have > twice the amount of coverage compared to genome")

    def gene_sample_filtered_on_quality(self):

        data = self.gene_sample_df

        breadth_field = "breadth_{depth}x".format(depth=self.threshold_depth)
        data = data[data[breadth_field] > self.threshold_breadth]

        # data = data[(data.AAcoverage_perc < 1.5)]
        # data = data[(data.AAcoverage_perc > 0.2)]

        return data

    def gene_sample_filtered_on_sample_coverage(self, data):
        # also filter out complete samples before showing pN/pS
        # we do not care for samples with a 1x horizontal coverage below 5%
        filtered_sample = self.sample_breadth_df[self.sample_breadth_df.breadth_1x.ge(0.05)]

        data = data.merge(filtered_sample,
                          left_on=data["sample"],
                          right_on=filtered_sample["sample"],
                          how="inner").drop(["key_0"], axis=1)

        data.rename(columns={'sample_x': 'sample'}, inplace=True)

        return data

    # https://seaborn.pydata.org/generated/seaborn.heatmap.html
    def create_heatmap(self, data, measure, title):

        data = data[['Protein', 'sample', measure]]
        data = data.sort_values("Protein")
        data = data.set_index(["sample", "Protein"])

        # nr_of_samples = 3
        # nr_of_genes = 10
        # data = data.tail(nr_of_samples * nr_of_genes)

        # convert from multi-index to cross-product table
        data = data.unstack()

        # data = data.transpose() # if you would like to switch columns and rows

        # rename columns, unstack
        data.columns = ["_".join(x) for x in data.columns.ravel()]

        # data.columns = [x.split("_")[-1] for x in data.columns]
        data.columns = ["_".join(x.split("_")[3:]) for x in data.columns]

        plt.clf()  # clear the current figure (always do this before any new plot)
        ax = sns.heatmap(data, cmap="seismic", annot=False)
        plt.title(title + " ({breadth}/{depth}x)".
                  format(depth=self.threshold_depth, breadth=self.threshold_breadth))

        figure_name = "{}{}gene_plots.heat_map.{breadth}.{depth}x.pdf".format(
            self.plot_dir, self.dir_sep, measure.replace("/", "_"),
            depth=self.threshold_depth, breadth=self.threshold_breadth
        )
        plt.savefig(figure_name)

    # we want to see what genes have the highest and lowest pN/pS scores
    # based on self.gene_sample_df
    # that can be further aggregated into self.gene_df
    # and then merged with self.gene_annotation.df for annotation
    def score_genes_and_output_ordered_genes_to_file(self):

        filtered_gene_sample_df = self.gene_sample_filtered_on_quality()

        self.gene_df = filtered_gene_sample_df.groupby("Protein").agg(
            {
                'log10_pN/pS': ["mean", "count", "std"]
            ,   'entropy_mean': ["mean", "count"]
            }).reset_index()

        self.gene_df.columns = ["_".join(x) for x in self.gene_df.columns.ravel()]
        self.gene_df.rename(columns={'Protein_': 'Protein'}, inplace=True)

        self.gene_df = self.gene_df.sort_values(by='log10_pN/pS_mean', ascending=False)

        merge_df = self.gene_df.merge(self.gene_anno_df
                                     , left_on=self.gene_df.Protein
                                     , right_on=self.gene_anno_df.Protein
                                     , how='inner')
        merge_df.rename(columns={'key_0': 'Protein'}, inplace=True)

        merge_df = merge_df[['Protein','log10_pN/pS_mean','log10_pN/pS_std','log10_pN/pS_count','Annotation']]

        filename = self.plot_dir + self.dir_sep + "crassphage_pN_pS_values.{breadth}.{depth}x.txt".format(
            depth=self.threshold_depth, breadth=self.threshold_breadth
        )
        merge_df.to_csv(filename, index=False, sep='\t')

    def score_samples(self):
        filtered_gene_sample_df = self.gene_sample_filtered_on_quality()

        self.sample_df = filtered_gene_sample_df.groupby("sample").agg(
            {
                'log10_pN/pS':      ["mean", "count", "std"]
            ,   'entropy_mean':     ["mean", "count"]
            ,   'AAcoverage_mean':  ["mean"]
            }).reset_index()

        self.sample_df.columns = ["_".join(x) for x in self.sample_df.columns.ravel()]
        self.sample_df.rename(columns={'sample_': 'sample'}, inplace=True)

        self.sample_df = self.sample_df.sort_values(by='entropy_mean_mean', ascending=False)

    def create_box_plots(self, min_nr_samples=3):

        # box plots for log10_pN/pS

        # filter out the genes that do not have a minimal presence of min_nr_of_samples occurences in the samples
        filter_data = self.gene_df[self.gene_df['log10_pN/pS_count'] >= min_nr_samples]

        # head only works because self.gene_df is already ordered by log10_pN/pS_mean in score_genes_..()
        top10_data = filter_data.head(10)[['Protein', 'log10_pN/pS_mean']]
        self.create_box_plot(top10_data, "Protein", "log10_pN/pS",
                             "top 10 pos selection present in at least {} samples".format(min_nr_samples))

        bottom10_data = filter_data.tail(10)[['Protein', 'log10_pN/pS_mean']]
        self.create_box_plot(bottom10_data, "Protein", "log10_pN/pS",
                             "top 10 most conserved present in at least {} samples".format(min_nr_samples))

        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "Protein", "log10_pN/pS",
                             "top and bottom 10 present in at least {} samples".format(min_nr_samples))

        self.create_box_plot(self.gene_df, "Protein", "log10_pN/pS", "all genes")

        # box plots for ENTROPY
        filter_data = self.gene_df[self.gene_df['entropy_mean_count'] >= min_nr_samples]

        top10_data = filter_data.head(10)[['Protein', 'entropy_mean_mean']]
        self.create_box_plot(top10_data, "Protein", "entropy_mean",
                             "top 10 internal var present in at least {} samples".format(min_nr_samples))

        bottom10_data = filter_data.tail(10)[['Protein', 'entropy_mean_mean']]
        self.create_box_plot(bottom10_data, "Protein", "entropy_mean",
                             "bottom 10 internal var present in at least {} samples".format(min_nr_samples))

        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "Protein", "entropy_mean",
                             "top and bottom 10 internal var in at least {} samples".format(min_nr_samples))

        # new aggregation per sample (min_nr_samples should now be read as min_nr_genes)
        filter_data = self.sample_df[self.sample_df['log10_pN/pS_count'] >= min_nr_samples]

        top10_data = filter_data.head(10)
        bottom10_data = filter_data.tail(10)
        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "sample", "entropy_mean",
                             "top and bottom 10 entropy with at least {} genes".format(min_nr_samples))

        filter_data = filter_data.sort_values(by='AAcoverage_mean_mean', ascending=False)

        top10_data = filter_data.head(10)
        self.create_box_plot(top10_data, "sample", "AAcoverage_mean",
                             "top 10 coverage with at least {} genes".format(min_nr_samples))

        bottom10_data = filter_data.tail(10)
        self.create_box_plot(bottom10_data, "sample", "AAcoverage_mean",
                             "bottom 10 coverage with at least {} genes".format(min_nr_samples))

        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "sample", "AAcoverage_mean",
                             "top and bottom 10 coverage with at least {} genes".format(min_nr_samples))

    def create_box_plot(self, filter_data, agg_field, measure, title):

        data = self.gene_sample_filtered_on_quality()
        # data = self.gene_sample_df

        data = data.merge(filter_data
                             , left_on=data[agg_field]
                             , right_on=filter_data[agg_field]
                             , how='inner')

        data.rename(columns={"{}_x".format(agg_field): agg_field}, inplace=True)

        data = data[[agg_field, measure]]

        # add mean of measure per gene and then merge with the original dataset
        # in order to sort the genes on the mean of the measure
        grouped = data.groupby(agg_field).mean()

        data = data.merge(grouped
                          , left_on=data[agg_field]
                          , right_on = grouped.index
                          , how='inner').sort_values(["{}_y".format(measure), agg_field], ascending=False).\
            rename(columns={"{}_x".format(measure): measure })

        plt.clf()
        plt.title(title + " ({breadth}/{depth}x)".
                  format(depth=self.threshold_depth, breadth=self.threshold_breadth))

        sns.set(style="ticks")

        # to do: We get a warning on the percentile calculations (implicit in box plot) for the infinite values
        # we should probably recalculate p_N/p_S with a pseudocount
        sns.boxplot(x=measure, y=agg_field, data=data,
                    whis="range", palette="vlag")

        sns.swarmplot(x=measure, y=agg_field, data=data,
                      size=2, color=".3", linewidth=0)

        figure_name = "{}{}gene_plots.{}.box_plot.{}.{}.{breadth}.{depth}x.pdf".format(
            self.plot_dir, self.dir_sep, agg_field, measure.replace("/", "_"), title.replace(" ", "_"),
            depth=self.threshold_depth, breadth=self.threshold_breadth
        )
        plt.savefig(figure_name)

    def create_line_plots_for_pn_ps(self):
        # make line plots based on
        # aggregration of bin_sample_df
        # to bin_df

        nr_lines = len(self.bin_sample_df)

        filtered_bin_data = \
            self.bin_sample_df[(self.bin_sample_df.AAcoverage_perc > 0.2) & (self.bin_sample_df.AAcoverage_cv < 0.2)]

        nr_filtered_lines = len(filtered_bin_data)

        # to do: do not take mean of log but instead log of mean
        data = filtered_bin_data.groupby(["Protein", "AAPosition"]).agg(
            {
                'log10_pN_pS_60': ["mean"]
            }).reset_index()
        data.columns = ["_".join(x) for x in data.columns.ravel()]
        data.rename(columns={'Protein_': 'Protein'}, inplace=True)
        data.rename(columns={'AAPosition_': 'AAPosition'}, inplace=True)

        data.rename(columns={'log10_pN_pS_60_mean': 'log10_pN_pS_60'}, inplace=True)

        data = data.join(data.groupby('Protein')['AAPosition'].max(), on='Protein', rsuffix='_max')

        self.line_plots_for_measure(data, "log10_pN_pS_60")

    def line_plots_for_measure(self, data, measure, ylim_bottom=None):

        genes = data.Protein.unique()

        # determine the number of genes and make a plot for every set of 6 genes
        nr_plots = int(np.round((len(genes)/6)))

        # we loop through the plots and subplots and then select the next gene instead of looping through the genes
        i = 0
        # we always make one extra plot (so some of the sub plots of the last plot may be empty)
        for j in range(0,nr_plots):

            fig, axs = plt.subplots(3,2)
            fig.tight_layout(pad=2)
            for row in axs:
                for ax in row:
                    #to prevent
                    if len(genes) > i:
                        gene = genes[i]
                        i = i + 1
                        gene_data = data[data.Protein == gene]

                        x_max = gene_data.AAPosition_max.max()

                        sns.lineplot(legend=None, data=gene_data, x="AAPosition", y=measure, ax=ax)
                        ax.set_ylim(bottom=ylim_bottom)
                        ax.set_xlim(left=1)
                        ax.set_xlim(right=x_max)
                        ax.set_title(gene)

            start = 1 + j*6
            end = 6 + j*6
            plot_name = self.plot_dir + self.dir_sep + "{}_{}_{}.png".format(measure, start, end)

            plt.savefig(plot_name)
            plt.close()

    def do_analysis(self, min_nr_samples, ref):

        self.read_files(ref)

        self.create_plot_dir(ref)

        self.create_line_plots_for_pn_ps()

        self.prepare_data_for_plots()

        self.create_histograms()

        self.create_heatmaps()

        self.score_genes_and_output_ordered_genes_to_file()

        self.score_samples()

        # create box plots for top 10 and bottom 10 pN/pS values for genes with a minimum nr of samples for that gene
        # after applying the quality filter
        self.create_box_plots(min_nr_samples)


def do_analysis(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)
    parser.add_argument("-r", "--ref", dest="ref",
                        help="reference genome id", metavar="[ref]", required=True)

    parser.add_argument("-rd", "--ref_dir", dest="ref_dir",
                        help="directory reference genomes", metavar="[ref_dir]", required=True)

    parser.add_argument("-ns", "--min_nr_samples", dest="min_nr_samples", metavar="[min_nr_samples]", type=int,
                        help="minimal nr of samples for genes")

    parser.add_argument("-td", "--threshold_depth", dest="threshold_depth", metavar="[threshold_depth]", type=int,
                        help="threshold for the depth, 1,10 or 50, will be combined with threshold for breadth")

    parser.add_argument("-tb", "--threshold_breadth", dest="threshold_breadth", metavar="[threshold_breadth]",
                        type=float,
                        help="threshold for breadth, between 0 and 1, will be combined with threshold for depth")

    args = parser.parse_args(args_in)

    print("Start running MakeGenePlots")
    print("sample_dir: " + args.sample_dir)
    print("reference gemome id (ref): " + args.ref)

    if not args.min_nr_samples:
        args.min_nr_samples = 3
    print("minimal nr of samples = {}".format(args.min_nr_samples))
    if not args.threshold_depth:
        args.threshold_depth = 10
    if not args.threshold_breadth:
        args.threshold_breadth = 0.95
    print("threshold depth   : {depth}".format(depth=args.threshold_depth))
    print("threshold breadth : {breadth}".format(breadth=args.threshold_breadth))

    make = MakeGenePlots(args.sample_dir, args.ref_dir, args.threshold_depth, args.threshold_breadth)

    make.do_analysis(args.min_nr_samples, args.ref)


if __name__ == "__main__":
    do_analysis(sys.argv[1:])

#TODO for testing, do not use in production
# sample_dir=r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"
# ref="crassphage_refseq"
# # ref="sib1_ms_5"
# rd = r"D:\17 Dutihl Lab\source\phages_scripts\mgx\ref_seqs"
# do_analysis(["-d", sample_dir, "-rd", rd, "-r", ref, "-ns", "1", "-td", "10", "-tb", "0.95"])


