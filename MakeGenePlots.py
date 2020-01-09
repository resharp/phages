import argparse
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns;sns.set()
import sys

#MakeGenePlots
# create multiple heatmaps
#   fst with the log pN/pS values
#   snd with the missing genes (just 0/1 for clarity)
#   trd distribution plots [integration over structural genes or other categories]
#
#   all sample_gene measures are in small separate files
#   we aggregate them to gene level (aggregation over all samples)
#   therefore we load all files and merge them into one dataframe
#
#masking missing values in heatmaps:
#https://github.com/mwaskom/seaborn/issues/375

class MakeGenePlots:

    sample_dir = ""
    dir_sep = ""
    plot_dir = ""

    gene_sample_df = None
    gene_df = None
    gene_anno_df = None

    def __init__(self, sample_dir):

        self.sample_dir = sample_dir
        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'MakeGenePlots.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.INFO)

    def read_files(self):

        logging.info("start reading tables")
        self.read_gene_sample_measures()

        self.read_gene_annotation()

        logging.info("finished reading tables")

    def read_gene_sample_measures(self):

        subfolders = [f.path for f in os.scandir(self.sample_dir) if f.is_dir()]
        samples = []

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            sample_name = subfolder + self.dir_sep + sample + "_gene_measures.txt"

            if os.path.isfile(sample_name) and os.path.getsize(sample_name) > 0:
                logging.info("processing {}".format(sample_name))

                gene_sample_df = pd.read_csv(sample_name
                                         ,  sep='\t'
                                         )
                samples.append(gene_sample_df)

        self.gene_sample_df = pd.concat(samples)

    def read_gene_annotation(self):

        #TODO: also make the annotation file a parameter of MakeGenePlots
        gene_anno_file_name = self.sample_dir + self.dir_sep + "crassphage_gene_list.txt"

        self.gene_anno_df = pd.read_csv(gene_anno_file_name
                                        ,   sep='\t'
                                        ,   header=None
                                        ,   usecols=[0,4]
                                        ,   names=["Protein", "Annotation"]
                                        )
        pass

    def prepare_data_for_plots(self):

        #here we can prepare the data further for displaying purposes
        #it could also have been done in CalcDiversiMeasures

        #we want to see the missing genes, and missing genes are if there are no mappings or hardly any
        self.gene_sample_df['missing_gene'] = 0
        self.gene_sample_df.loc[self.gene_sample_df.AAcoverage_perc < 0.02 , 'missing_gene'] = 1
        self.gene_sample_df.loc[np.isnan(self.gene_sample_df.AAcoverage_perc), 'missing_gene'] = 1

        self.gene_sample_df['double_coverage'] = 0
        self.gene_sample_df.loc[self.gene_sample_df.AAcoverage_perc > 2 , 'double_coverage'] = 1

        # data_debug = self.gene_df[['sample','Protein', 'AAcoverage_perc', 'missing_gene', 'double_coverage']]

        #make binary plot for positive selection or conservation
        self.gene_sample_df['positive_selection'] = 0
        self.gene_sample_df.loc[self.gene_sample_df["log10_pN/pS"] > 0.1, 'positive_selection'] = 1
        self.gene_sample_df.loc[self.gene_sample_df["log10_pN/pS"] < -0.1, 'positive_selection'] = -1

    def create_plot_dir(self):

        self.plot_dir = self.sample_dir + self.dir_sep + "GenePlots"
        os.makedirs(self.plot_dir, exist_ok=True)

    def create_histograms(self):

        #create a histogram of the occurence of genes across samples
        data = self.gene_sample_filtered_on_quality()

        #we only need the counts per gene

        counts_df = data.groupby("Protein").agg(
            {
                'log10_pN/pS': ["count"]
            }).reset_index()
        counts_df.columns = ["_".join(x) for x in counts_df.columns.ravel()]
        counts_df.rename(columns={'Protein_': 'Protein'}, inplace=True)

        min = counts_df['log10_pN/pS_count'].min()
        max = counts_df['log10_pN/pS_count'].max()

        sns.distplot(counts_df[['log10_pN/pS_count']], bins=(max-min), kde=False)

        plt.title("Distribution of sample presence for genes")
        figure_name = "{}{}gene_plots.gene_sample_counts.pdf".format(self.plot_dir, self.dir_sep)

        plt.savefig(figure_name)

    def create_heatmaps(self):

        #filter on quality
        data = self.gene_sample_filtered_on_quality()

        # make a heatmap of the log10_pN/pS based on multiple samples
        self.create_heatmap(data, "log10_pN/pS", "Log 10 of pN/pS (blue = positive selection)")
        self.create_heatmap(data, "positive_selection", "log10_pN/pS either > 0.1 or < - 0.1")

        self.create_heatmap(data, "SndAAcnt_perc_polymorphism_mean", "Within sample AA variation in genes")

        self.create_heatmap(data, "AAcoverage_perc", "Coverage percentage (compared to whole genome)")

        self.create_heatmap(data, "entropy_mean", "Mean codon entropy (base 10)")

        #make a heatmap of quality measure AAcoverage_cv
        self.create_heatmap(data, "AAcoverage_cv", "Internal coefficient of variation per gene")

        #use unfiltered data for quality measures
        data = self.gene_sample_df
        self.create_heatmap(data, "missing_gene", "Missing genes (based on coverage < 2% rest genome)")

        self.create_heatmap(data, "double_coverage", "Genes that have > twice the amount of coverage compared to genome")

    def gene_sample_filtered_on_quality(self):

        data = self.gene_sample_df
        data = data[data.AAcoverage_cv < 0.2]

        # data = data[(data.AAcoverage_perc < 1.5)]
        data = data[(data.AAcoverage_perc > 0.2)]

        #TODO: should we also exclude samples with few values for genes left after filtering?
        #first make distribution plot of # of genes per sample

        return data

    #https://seaborn.pydata.org/generated/seaborn.heatmap.html
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
        ax = sns.heatmap(data, cmap="YlGnBu", annot=False)
        plt.title(title)

        figure_name = "{}{}gene_plots.heat_map.{}.pdf".format(self.plot_dir, self.dir_sep, measure.replace("/", "_"))

        plt.savefig(figure_name)

    #we want to see what genes have the highest and lowest pN/pS scores
    #based on self.gene_sample_df
    #that can be further aggregated into self.gene_df
    #and then merged with self.gene_annotation.df for annotation
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

        merge_df.to_csv(self.plot_dir + self.dir_sep + "crassphage_pN_pS_values.txt", index=False, sep='\t')

    def create_box_plots(self, min_nr_samples=3):

        #box plots for log10_pN/pS

        #filter out the genes that do not have a minimal presence of min_nr_of_samples occurences in the samples
        filter_data = self.gene_df[self.gene_df['log10_pN/pS_count'] >= min_nr_samples]

        top10_data = filter_data.head(10)[['Protein', 'log10_pN/pS_mean']]
        self.create_box_plot(top10_data, "log10_pN/pS",
                             "top 10 pos selection present in at least {} samples".format(min_nr_samples))

        bottom10_data = filter_data.tail(10)[['Protein', 'log10_pN/pS_mean']]
        self.create_box_plot(bottom10_data, "log10_pN/pS",
                             "top 10 most conserved present in at least {} samples".format(min_nr_samples))

        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "log10_pN/pS",
                             "top and bottom 10 present in at least {} samples".format(min_nr_samples))

        self.create_box_plot(self.gene_df, "log10_pN/pS", "all genes")

        #box plots for ENTROPY
        filter_data = self.gene_df[self.gene_df['entropy_mean_count'] >= min_nr_samples]

        top10_data = filter_data.head(10)[['Protein', 'entropy_mean_mean']]
        self.create_box_plot(top10_data, "entropy_mean",
                             "top 10 internal var present in at least {} samples".format(min_nr_samples))

        bottom10_data = filter_data.tail(10)[['Protein', 'entropy_mean_mean']]
        self.create_box_plot(bottom10_data, "entropy_mean",
                             "bottom 10 internal var present in at least {} samples".format(min_nr_samples))

        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "entropy_mean",
                             "top and bottom 10 internal var in at least {} samples".format(min_nr_samples))


    def create_box_plot(self, filter_data, measure, title):

        data = self.gene_sample_filtered_on_quality()
        # data = self.gene_sample_df

        data = data.merge(filter_data
                             , left_on=data.Protein
                             , right_on=filter_data.Protein
                             , how='inner')

        data.rename(columns={"Protein_x": "Protein"}, inplace=True)

        data = data[["Protein", measure]]

        #add mean of measure per gene and then merge with the original dataset
        # in order to sort the genes on the mean of the measure
        grouped = data.groupby("Protein").mean()

        data = data.merge(grouped
                          , left_on=data.Protein
                          , right_on = grouped.index
                          , how='inner').sort_values(["{}_y".format(measure), "Protein"], ascending=False).\
            rename(columns={"{}_x".format(measure): measure })

        plt.clf()
        plt.title(title)
        sns.set(style="ticks")

        #TODO: We get a warning on the percentile calculations (implicit in box plot) for the infinite values
        #we should probably recalculate p_N/p_S with a pseudocount
        sns.boxplot(x=measure, y="Protein", data=data,
                    whis="range", palette="vlag")

        sns.swarmplot(x=measure, y="Protein", data=data,
                      size=2, color=".3", linewidth=0)

        figure_name = "{}{}gene_plots.box_plot.{}.{}.pdf".format(self.plot_dir, self.dir_sep,
                                                                 measure.replace("/", "_"),
                                                                 title.replace(" ", "_"))
        # plt.show()
        plt.savefig(figure_name)

    def do_analysis(self, min_nr_samples):

        self.read_files()

        self.prepare_data_for_plots()

        self.create_plot_dir()

        self.create_histograms()

        self.create_heatmaps()

        self.score_genes_and_output_ordered_genes_to_file()

        #TODO in same way: what samples have the highest entropy?

        #create box plots for top 10 and bottom 10 pN/pS values for genes with a minimum nr of samples for that gene
        #after applying the quality filter
        self.create_box_plots(min_nr_samples)

def do_analysis(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    parser.add_argument("-ns", "--min_nr_samples", dest="min_nr_samples", metavar="[min_nr_samples]", type=int,
                        help="minimal nr of samples for genes")

    args = parser.parse_args(args_in)

    print("Start running MakeGenePlots")
    print("sample_dir: " + args.sample_dir)

    if not args.min_nr_samples:
        args.min_nr_samples = 3
    print("minimal nr of samples = {}".format(args.min_nr_samples))

    make = MakeGenePlots(args.sample_dir)

    make.do_analysis(args.min_nr_samples)

if __name__ == "__main__":
    do_analysis(sys.argv[1:])

#TODO for testing, do not use in production
# sample_dir = r"D:\17 Dutihl Lab\_tools\diversitools"
# do_analysis(["-d", sample_dir, "-ns", "4"])


