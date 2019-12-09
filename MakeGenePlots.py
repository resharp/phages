import logging
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns;sns.set()

if os.name == "nt":
    sample_dir = r"D:\17 Dutihl Lab\_tools\diversitools"
else:
    sample_dir = "/hosts/linuxhome/mutant31/tmp/richard/crassphage_samples"

#MakeGenePlots
# create two heatmaps
#   fst with the log dN/ds values
#   snd with the missing genes (just 0/1 for clarity? or use the normalized coverage?)
#   trd distribution plots [integration over structural genes or other categories]
#
#   all sample_gene measures are in small separate files
#   we will aggregate them to gene level (aggregation over all samples)
#   therefore we load all files and merge them into one dataframe

#masking missing values in heatmaps
#https://github.com/mwaskom/seaborn/issues/375

class MakeGenePlots:

    sample_dir = ""
    dir_sep = ""

    gene_df = None
    # gene_sample_df = None.
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

    def make_heatmaps(self):

        #https://seaborn.pydata.org/generated/seaborn.heatmap.html

        # make a heatmap of the log10_dN/dS based on two (or more) samples
        # first load multiple samples in one file
        # we might also use cat for that?

        data = self.gene_df[['Protein', 'sample', 'log10_dN/dS']]

        data = data.sort_values("Protein")

        data = data.set_index(["sample", "Protein"])

        # data = data.head(40)

        data = data.unstack()

        # data = data.transpose() # if you would like to switch columns and rows

        plt.clf() #clear the current figure (always do this before any new plot)
        ax = sns.heatmap(data,  cmap="YlGnBu", annot=False)
        # plt.show()

        figure_name = "{}{}gene_plots.{}.pdf" \
            .format(self.sample_dir, self.dir_sep, "log10_dN_dS")

        plt.show()
        #plt.savefig(figure_name)


make = MakeGenePlots(sample_dir)

make.read_files()
make.make_heatmaps()

