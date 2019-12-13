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

    gene_df = None

    gene_categories_df = None

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

        self.gene_df = pd.concat(samples)

        logging.debug("finished reading tables")


    def prepare_data(self):

        #here we can prepare the data further for displaying purposes
        #it could also have been done in CalcDiversiMeasures

        #we want to see the missing genes, and missing genes are if there are no mappings or hardly any
        self.gene_df['missing_gene'] = 0
        self.gene_df.loc[self.gene_df.AAcoverage_perc < 0.02 , 'missing_gene'] = 1
        self.gene_df.loc[np.isnan(self.gene_df.AAcoverage_perc), 'missing_gene'] = 1

        self.gene_df['double_coverage'] = 0
        self.gene_df.loc[self.gene_df.AAcoverage_perc > 2 , 'double_coverage'] = 1

        # data_debug = self.gene_df[['sample','Protein', 'AAcoverage_perc', 'missing_gene', 'double_coverage']]

        #make binary plot for positive selection or conservation
        self.gene_df['positive_selection'] = 0
        self.gene_df.loc[self.gene_df["log10_pN/pS"] > 0.1, 'positive_selection'] = 1
        self.gene_df.loc[self.gene_df["log10_pN/pS"] < -0.1, 'positive_selection'] = -1

    def create_plot_dir(self):

        self.plot_dir = self.sample_dir + self.dir_sep + "GenePlots"
        os.makedirs(self.plot_dir, exist_ok=True)

    #https://seaborn.pydata.org/generated/seaborn.heatmap.html
    def create_heatmaps(self):

        self.create_plot_dir()

        #filter on quality
        data = self.filter_on_quality(self.gene_df)

        # make a heatmap of the log10_dN/dS based on multiple samples
        self.create_heatmap(data, "log10_pN/pS", "Log 10 of pN/pS (blue = positive selection)")
        self.create_heatmap(data, "positive_selection", "log10_pN/pS either > 0.1 or < - 0.1")

        self.create_heatmap(data, "SndAAcnt_perc_filtered_mean", "Within sample AA variation in genes")

        self.create_heatmap(data, "AAcoverage_perc", "Coverage percentage (compared to whole genome)")

        #make a heatmap of quality measure AAcoverage_cv
        self.create_heatmap(data, "AAcoverage_cv", "Internal coefficient of variation per gene")

        #use unfiltered data for quality measures
        data = self.gene_df
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

def make_plots(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    args = parser.parse_args(args_in)

    print("Start running MakeGenePlots")
    print("sample_dir: " + args.sample_dir)

    make = MakeGenePlots(args.sample_dir)

    make.read_files()

    make.prepare_data()

    make.create_heatmaps()

if __name__ == "__main__":
    make_plots(sys.argv[1:])

#TODO for testing, do not use in production
# sample_dir = r"D:\17 Dutihl Lab\_tools\diversitools"
# make_plots(["-d", sample_dir])
