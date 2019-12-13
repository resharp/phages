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
#   fst with the log dN/ds values
#   snd with the missing genes (just 0/1 for clarity? or use the normalized coverage?)
#   trd distribution plots [integration over structural genes or other categories]
#
#   all sample_gene measures are in small separate files
#   we aggregate them to gene level (aggregation over all samples)
#   therefore we load all files and merge them into one dataframe

#masking missing values in heatmaps
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
                            level=logging.DEBUG)

    def read_files(self):

        logging.debug("start reading tables")
        self.read_gene_sample_measures()

        #TODO: add read gene_annotaion
        self.read_gene_annotation()

        logging.debug("finished reading tables")

    def read_gene_sample_measures(self):

        subfolders = [f.path for f in os.scandir(self.sample_dir) if f.is_dir()]
        samples = []

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            sample_name = subfolder + self.dir_sep + sample + "_gene_measures.txt"

            if os.path.isfile(sample_name) and os.path.getsize(sample_name) > 0:
                logging.debug("processing {}".format(sample_name))

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

    #https://seaborn.pydata.org/generated/seaborn.heatmap.html
    def create_heatmaps(self):

        #filter on quality
        data = self.filter_on_quality(self.gene_sample_df)

        # make a heatmap of the log10_dN/dS based on multiple samples
        self.create_heatmap(data, "log10_pN/pS", "Log 10 of pN/pS (blue = positive selection)")
        self.create_heatmap(data, "positive_selection", "log10_pN/pS either > 0.1 or < - 0.1")

        self.create_heatmap(data, "SndAAcnt_perc_filtered_mean", "Within sample AA variation in genes")

        self.create_heatmap(data, "AAcoverage_perc", "Coverage percentage (compared to whole genome)")

        #make a heatmap of quality measure AAcoverage_cv
        self.create_heatmap(data, "AAcoverage_cv", "Internal coefficient of variation per gene")

        #use unfiltered data for quality measures
        data = self.gene_sample_df
        self.create_heatmap(data, "missing_gene", "Missing genes (based on coverage < 2% rest genome)")

        self.create_heatmap(data, "double_coverage", "Genes that have > twice the amount of coverage compared to genome")

    def filter_on_quality(self, data):

        data = data[data.AAcoverage_cv < 0.2]

        data = data[(data.AAcoverage_perc < 1.5)]
        data = data[(data.AAcoverage_perc > 0.5)]
        return data

    def create_heatmap(self, data, measure_field, title):

        data = data[['Protein', 'sample', measure_field]]
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

        # can we remove gp02 and gp03?
        # data = data.drop(["gp02", "gp03", "gp04"], axis=1)

        plt.clf()  # clear the current figure (always do this before any new plot)
        ax = sns.heatmap(data, cmap="YlGnBu", annot=False)
        plt.title(title)

        figure_name = "{}{}gene_plots.{}.pdf".format(self.plot_dir, self.dir_sep, measure_field.replace("/", "_"))

        #plt.show()
        plt.savefig(figure_name)

    #now we want to see what genes have the highest and lowest pN/pS scores
    #based on self.gene_sample_df
    #that can be further aggregated into self.gene_df
    #and then merged with self.gene_annotation.df for annotation
    def score_genes(self):

        self.gene_df = self.gene_sample_df.groupby("Protein").agg(
            {
                'log10_pN/pS': ["mean", "count", "std"]
            }).reset_index()

        self.gene_df.columns = ["_".join(x) for x in self.gene_df.columns.ravel()]
        self.gene_df.rename(columns={'Protein_': 'Protein'}, inplace=True)

        self.gene_df = self.gene_df.sort_values(by='log10_pN/pS_mean', ascending=False)

        merge_df = self.gene_df.merge(self.gene_anno_df
                                     , left_on=self.gene_df.Protein
                                     , right_on=self.gene_anno_df.Protein
                                     , how='inner')
        merge_df.rename(columns={'key_0': 'Protein'}, inplace=True)

        merge_df = merge_df[['Protein','log10_pN/pS_mean','log10_pN/pS_std','Annotation']]

        merge_df.to_csv(self.plot_dir + self.dir_sep + "crassphage_pN_pS_values.txt", index=False, sep='\t')

    def do_analysis(self):

        self.read_files()

        self.prepare_data_for_plots()

        self.create_plot_dir()

        self.create_heatmaps()

        self.score_genes()

def do_analysis(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    args = parser.parse_args(args_in)

    print("Start running MakeGenePlots")
    print("sample_dir: " + args.sample_dir)

    make = MakeGenePlots(args.sample_dir)

    make.do_analysis()

if __name__ == "__main__":
    do_analysis(sys.argv[1:])

#TODO for testing, do not use in production
# sample_dir = r"D:\17 Dutihl Lab\_tools\diversitools"
# do_analysis(["-d", sample_dir])
