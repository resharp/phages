import argparse
import logging
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns;sns.set()


# we will calculate per gene
#   average AA coverage
#   standard dev
#   normalized stdev (coefficient of variation)
#   sum(CntNonSys)
#   sum(CntSyn)
#   pN/pS
#   % TopCodoncnt
#   % SndCodoncnt
#   % TrdCodoncnt
#   entropy
#   ..
#
#   we do this on one sample and output a file with measures per gene
#   there is no filtering here
class CalcDiversiMeasures:

    aa_table_name = ""
    codon_table_name = ""
    gene_table_name = ""

    sample_dir = ""
    plot_dir = ""

    dir_sep = ""
    aa_df = None
    gene_df = None
    codon_df = None

    def __init__(self, sample_dir):

        logging.basicConfig(filename='CalcDiversiMeasures.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.DEBUG)

        self.sample_dir = sample_dir

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

    def read_files(self, sample):
        logging.debug("start reading tables")

        self.read_aa_table(sample)

        self.read_codon_table()

        self.integrate_tables()

        logging.debug("finished reading tables")

    def read_aa_table(self, sample):

        self.aa_table_name = self.sample_dir + self.dir_sep + sample + self.dir_sep + sample + "_AA_clean.txt"
        self.aa_df = pd.read_csv(self.aa_table_name
                                 ,  sep='\t'
                                 ,  usecols=range(2,26)
                                 )
        # automatic columns with first two skipped
        # SKIPPED: Sample	Chr
        # Protein	AAPosition	RefAA	RefSite	RefCodon	FstCodonPos	SndCodonPos	TrdCodonPos
        # CntNonSyn	CntSyn	NbStop	TopAA	TopAAcnt	SndAA	SndAAcnt	TrdAA	TrdAAcnt
        # TopCodon	TopCodoncnt	SndCodon	SndCodoncnt	TrdCodon	TrdCodoncnt	AAcoverage

    def read_codon_table(self):
        self.codon_table_name = self.sample_dir + self.dir_sep + "codon_syn_non_syn_probabilities.txt"

        self.codon_df = pd.read_csv(self.codon_table_name
                                    ,   sep=","
                                    ,   usecols=[0,2,3,4])

    def integrate_tables(self):

        # also interesting: check codon usage signature in genes?
        # because we do not want our measures be biased for codon usage

        self.aa_df = self.aa_df.merge(self.codon_df
                                    , left_on=self.aa_df.RefCodon
                                     , right_on=self.codon_df.codon
                                     , how='left')#.reset_index()

        #todo: hoe kan refcodon nan zijn?
        # merge_df_short = merge_df[['Protein','AAPosition','RefCodon', 'codon', 'syn', 'non_syn']]

    def calc_measures(self):

        # ENTROPY determine estimation of entropy per AA locus, based on three codon varieties
        self.aa_df['TopAAcnt_perc'] = (self.aa_df["TopAAcnt"] / self.aa_df["AAcoverage"])
        self.aa_df["SndAAcnt_perc"] = (self.aa_df["SndAAcnt"] / self.aa_df["AAcoverage"])
        self.aa_df['TrdAAcnt_perc'] = (self.aa_df["TrdAAcnt"] / self.aa_df["AAcoverage"])

        # we only keep the nans if the entropy for TopAAcnt_perc would be nan (=missing depth coverage)
        self.aa_df["entropy"] = np.abs(-self.aa_df['TopAAcnt_perc'] * np.log10(self.aa_df['TopAAcnt_perc'])
                                       -(self.aa_df['SndAAcnt_perc'] * np.log10(self.aa_df['SndAAcnt_perc'])).replace(np.nan,0)
                                       -(self.aa_df['TrdAAcnt_perc'] * np.log10(self.aa_df['TrdAAcnt_perc'])).replace(np.nan,0)
                                       )

        #DEFINITION of POLYMORPHISM
        #before aggregating first make a filtered derived measure for SndAAcnt_perc
        #we only want to keep percentages > 1% in the sample
        # SndAAcnt_perc_polymorphism
        self.aa_df['SndAAcnt_perc_polymorphism'] = 0
        self.aa_df.loc[self.aa_df["SndAAcnt_perc"] > 0.01, 'SndAAcnt_perc_polymorphism'] = self.aa_df["SndAAcnt_perc"]

        #total number of SNPs
        self.aa_df["CntSnp"] = self.aa_df["CntSyn"].replace(np.nan, 0) + self.aa_df["CntNonSyn"].replace(np.nan, 0)

        # debug_data = self.aa_df[["Protein", "AAPosition", "SndAAcnt_perc", "SndAAcnt_perc_polymorphism"]]

    def create_plot_dir(self, sample):

        self.plot_dir = self.sample_dir + self.dir_sep + sample + self.dir_sep + "DiversiMeasures"
        os.makedirs(self.plot_dir, exist_ok=True)

    def line_plots_for_coverage(self, sample):

        data = self.aa_df[["Protein", "AAPosition", "AAcoverage"]]
        genes = data.Protein.unique()

        # TODO: determine the number of genes and make a plot for every set of 6 genes
        nr_plots = int(np.round((len(genes)/6)))

        measure = "AAcoverage"

        i = 0
        for j in range(0,nr_plots):

            fig, axs = plt.subplots(3,2)
            fig.tight_layout(pad=2)
            for row in axs:
                for ax in row:
                    if len(genes) > i:
                        gene = genes[i]
                        i = i + 1
                        gene_data = data[data.Protein == gene]

                        sns.lineplot(legend=None, data=gene_data, x="AAPosition", y=measure, ax=ax)
                        ax.set_ylim(bottom=0)
                        ax.set_title(gene)

            start = 1 + j*6
            end = 6 + j*6
            plot_name = self.plot_dir + self.dir_sep + "coverage_genes_{}_{}.png".format(start, end)

            # plt.show()
            plt.savefig(plot_name)
            plt.close()

    def joint_plot_for_aa_changes_and_entropy(self):
        #here we show the AAPosition the x-axis and the entropy on the y-axis for every SCP
        #though entropy seems to be an attribute of a sample we would like to see if there is

        self.aa_df["SCP"] = 0
        self.aa_df.loc[self.aa_df.SndAAcnt_perc_polymorphism > 0.01, "SCP"] = 1

        data = self.aa_df[["Protein", "AAPosition", "TopCodon", "RefCodon", "SndAAcnt_perc_polymorphism", "SCP", "entropy"]]

        genes = data.Protein.unique()
        # for gene in genes:

        for gene in ["KP06_gp44", "KP06_gp45", "KP06_gp64"]:
            gene_data = data[data.Protein == gene]
            len_gene = len(gene_data)

            gene_data = gene_data[gene_data.SCP == 1]
            if len(gene_data) > 0:
                max_entropy = gene_data.entropy.max() + 0.02

                sns.set(style="darkgrid")
                tips = sns.load_dataset("tips")
                g = sns.jointplot("AAPosition", "entropy", data=gene_data,
                                  # kind="kde",
                                  xlim=(0, len_gene),
                                  ylim=(-0.01, np.where(max_entropy > 0, max_entropy, 0.5)),
                                  marginal_kws=dict(bins=np.round(len_gene/20).astype(int), rug=True),
                                  color="m", height=7)
                plt.title(gene)
                plt.show()

    def find_and_write_local_regions(self, sample):

        #now we want to determine the distances between SNPs in self.aa_df
        #and then take the 5 percentile shortest distances to identify regions that are close
        #for now we determine SNP as non-empty RefCodon unequal to non-empty TopCodon

        # SCP "single codon polymorphism" = difference between TopCodon and RefCodon
        scp_df = self.aa_df
        scp_df = scp_df[scp_df.TopCodon.notnull()]
        scp_df = scp_df[scp_df.RefCodon.notnull()]
        scp_df = scp_df[scp_df.TopCodon != scp_df.RefCodon]

        self.find_and_write_hypervariable_regions(sample, scp_df, "non syn regions compared to ref")

        #TODO: reconsider, finding high peaks of entropy might be interesting
        # however, it will not find broad ranges of high entropy (we have to find some balance for the height_percentile
        # and the perceptile in find_and_write_hypervariable_regions depending on our question)
        height_percentile = 95
        entropy_cutoff = np.percentile(self.aa_df[self.aa_df["entropy"].notnull()].entropy, height_percentile)

        entropy_df = self.aa_df[["Protein", "AAPosition", "SndAAcnt_perc_polymorphism", "entropy"]]
        entropy_df = entropy_df[(entropy_df.entropy > entropy_cutoff)]

        self.find_and_write_hypervariable_regions(sample, entropy_df, "high entropy regions")

        #we should also be able to do the same for any other positions
        #for example for distances between within-sample variations

        # scp_df = self.aa_df
        # scp_df = scp_df[scp_df.SndAAcnt_perc_polymorphism > 0.01]
        #
        # self.find_and_write_hypervariable_regions(sample, scp_df, "intra sample high diversity regions ")

    # poi_df contains Protein and positions of interest of filtered rows
    # poi_df only only needs to contain Protein and AAPosition and you can put anything in AAPosition
    # however, for debugging purposes, you might want to add filter fields like RefCodon and TopCodon
    # that have contributed to the determination of the positions of interest
    def find_and_write_hypervariable_regions(self, sample, poi_df, title):

        # poi_df = poi_df[["Protein", "AAPosition", "RefCodon", "TopCodon"]].reset_index()
        poi_df = poi_df[["Protein", "AAPosition"]].reset_index()

        #add some extra fields to use for positioning
        poi_df.rename(columns={'index': 'AAPositionGenome'}, inplace=True)
        poi_df['AAPositionGenome'] = poi_df['AAPositionGenome'] + 1 #+1 to correct for 0-based result of reset_index()

        poi_df = poi_df.reset_index()
        poi_df.rename(columns={'index': 'PoiPosition'}, inplace=True)
        poi_df['PoiPosition'] = poi_df['PoiPosition'] + 1 #+1 to correct for 0-based result of reset_index()

        #distance_bw is the backward distance
        poi_df['distance_bw'] = poi_df['AAPositionGenome'].diff()
        #distance_fw is the forward distance
        poi_df['distance_fw'] = poi_df['distance_bw'].shift(-1)

        #these next two lines are the core of the method: determine the 5% percentile
        #of the distance distribution between positions of interest
        percentile = 5
        distance_cutoff = np.percentile(poi_df[poi_df['distance_bw'].notnull()].distance_bw, percentile)

        #for experimenting with manually set distance cut-offs, not determined by percentile
        #TODO ************************
        #distance_cutoff = 1

        #apply filter on distance cutoff between POIs (positions-of-interest)
        poi_df = poi_df[(poi_df.distance_bw <= distance_cutoff) | (poi_df.distance_fw <= distance_cutoff)]

        #PoiDistance = the distance between FILTERED rows
        # => PoiDistance 1 means that these SCPs should be joined in a hypervariable region according to distance cut-off
        # e.g. when the cut-off is set at 2 the "distance" between rows can be 2, and still the PoiDistance is 1
        #first calculate the backward distance
        poi_df['PoiDistance_bw'] = poi_df['PoiPosition'].diff()
        #forward distance is just the backward distance of the next row
        poi_df['PoiDistance_fw'] = poi_df['PoiDistance_bw'].shift(-1)

        #now filter out just the joining regions
        poi_df = poi_df[(poi_df['PoiDistance_fw'] <= 1) | (poi_df['PoiDistance_bw'] <= 1)]

        #TODO: poi_df might also be used to calculate overlap over multiple samples

        gene_scp_df = poi_df.groupby('Protein')['AAPosition'].apply(list).reset_index(name='positions')
        # -> results in example of positions: [2,3,4,5,19,20], with multiple

        #for determining the joining regions, we start looping, because doing this with a set operation seems difficult
        data_for_df = []
        for index, row in gene_scp_df.iterrows():

            # example of positions: [2,3,4,5,19,20] -> should result in two regions [2,3,4,5] and [19,20]
            # if the distance_cutoff is < the distance between pos 5 and 19
            positions = row["positions"]
            #initialize previous position in such a way that the first condition will hold and a region will be started
            old_pos = positions[0] - distance_cutoff
            #always start a new region for each gene
            region = []
            for pos in positions:
                if (pos <= old_pos + distance_cutoff):
                    region.append(pos)
                else:
                    #first finish the region
                    if len(region) > 0:
                        first_aa_pos = region[0]
                        last_aa_pos = region[-1]
                        data_for_df.append([row["Protein"], first_aa_pos, last_aa_pos])
                    #then start a new region
                    region = []
                    region.append(pos)
                #save previous position
                old_pos = pos
            #finish the last region
            if len(region) > 0:
                first_aa_pos = region[0]
                last_aa_pos = region[-1]
                data_for_df.append([row["Protein"], first_aa_pos, last_aa_pos])

        columns = ["Protein", "AAPositionStart", "AAPositionEnd"]

        new_gene_poi_df = pd.DataFrame(columns=columns, data=data_for_df)

        new_gene_poi_df["length"] = new_gene_poi_df["AAPositionEnd"] - new_gene_poi_df["AAPositionStart"] + 1

        gene_region_table_name = self.sample_dir + self.dir_sep + sample + self.dir_sep + \
                                 sample + "_gene_" + title.replace(" ", "_") + ".txt"
        new_gene_poi_df.to_csv(gene_region_table_name, sep='\t', index=False)

    #aggregate measures
    def aggregate_measures(self, sample):

        #mean coverage of all genes

        #TODO add CntSnp median over all genes
        genome_coverage_mean = self.aa_df.AAcoverage.mean().round(decimals=4)

        self.gene_df = self.aa_df.groupby("Protein").agg(
            {
                'AAcoverage': ["mean", "std", "sum"],
                'TopAAcnt': ["mean", "std", "sum"],
                'SndAAcnt': ["mean", "std", "sum"],
                'TrdAAcnt': ["mean", "std", "sum"],
                'CntNonSyn': 'sum',
                'CntSyn': 'sum',
                'SndAAcnt_perc': 'mean',
                'SndAAcnt_perc_polymorphism': 'mean',
                'syn': 'sum',
                'non_syn': 'sum',
                'CntSnp': 'mean',
                'entropy': 'mean'
            }
        )

        #remove multi-index set on column axis
        #https://www.shanelynn.ie/summarising-aggregation-and-grouping-data-in-python-pandas/
        self.gene_df.columns = ["_".join(x) for x in self.gene_df.columns.ravel()]

        self.gene_df.AAcoverage_sum = pd.to_numeric(self.gene_df.AAcoverage_sum, downcast='unsigned', errors='coerce')
        self.gene_df.CntNonSyn_sum = pd.to_numeric(self.gene_df.CntNonSyn_sum, downcast='unsigned', errors='coerce')
        self.gene_df.CntSyn_sum = pd.to_numeric(self.gene_df.CntSyn_sum, downcast='unsigned', errors='coerce')
        self.gene_df.TopAAcnt_sum = pd.to_numeric(self.gene_df.TopAAcnt_sum, downcast='unsigned', errors='coerce')
        self.gene_df.SndAAcnt_sum = pd.to_numeric(self.gene_df.SndAAcnt_sum, downcast='unsigned', errors='coerce')
        self.gene_df.TrdAAcnt_sum = pd.to_numeric(self.gene_df.TrdAAcnt_sum, downcast='unsigned', errors='coerce')

        #derived measures
        self.gene_df["syn_ratio"] = self.gene_df["syn_sum"] / self.gene_df["non_syn_sum"]

        #take median of non-zero values over all amino acids (alle genes)
        count_snp_median = self.aa_df.CntSnp[self.aa_df["CntSnp"] != 0].median().round(decimals=4)
        snp_pseudo_count = np.sqrt(count_snp_median) / 2

        self.gene_df["pN/pS"] = self.gene_df["syn_ratio"] * \
            (self.gene_df["CntNonSyn_sum"] + snp_pseudo_count)/(self.gene_df["CntSyn_sum"] + snp_pseudo_count)

        self.gene_df["log10_pN/pS"] = np.where(self.gene_df["pN/pS"] > 0, np.log10(self.gene_df["pN/pS"]), 0)

        #quality measures
        self.gene_df["AAcoverage_perc"] = self.gene_df.AAcoverage_mean / genome_coverage_mean
        #coefficient of variation for a gene
        self.gene_df["AAcoverage_cv"] = self.gene_df.AAcoverage_std / self.gene_df.AAcoverage_mean

        #round all floats to four decimals
        self.gene_df = self.gene_df.round(decimals=4)

        #gene output table should contain sample (for later integrating over multiple samples)
        self.gene_df["sample"] = sample

    def write_aggregated_measures(self, sample):

        self.gene_table_name = self.sample_dir + self.dir_sep + sample + self.dir_sep + sample + "_gene_measures.txt"
        #the gene name is in the index
        self.gene_df.to_csv(path_or_buf=self.gene_table_name, sep='\t')

    def run_calc(self, sample):

        self.read_files(sample)

        self.calc_measures()

        self.create_plot_dir(sample)

        ##self.joint_plot_for_aa_changes_and_entropy()
        self.line_plots_for_coverage(sample)

        #TODO: only do this in "verbose" mode?
        self.find_and_write_local_regions(sample)

        self.aggregate_measures(sample)

        self.write_aggregated_measures(sample)

def run_calc(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    parser.add_argument("-s", "--sample", dest="sample",
                        help="sample with Diversitools input", metavar="[sample}", required=False)

    parser.add_argument("-a", "--all", action="store_true", default=False,
                        help="run all samples in sample directory", required=False)

    args = parser.parse_args(args_in)

    print("Start running CalcDiversiMeasures")
    print("sample_dir: " + args.sample_dir)
    if args.sample:
        print("sample: " + args.sample)
    if args.all:
        print("all samples: " + str(args.all))

    if (not args.all) and (not args.sample):
        print("Error: either explicitly specify -a for all samples or specify a specific sample with -s [sample]")
        return

    # some wrapper code that
    # determines if a sample name is given
    # or it is specified that all samples should be calculated

    calc = CalcDiversiMeasures(args.sample_dir)

    if args.sample:
        calc.run_calc(args.sample)
        return

    if args.all:
        subfolders = [f.path for f in os.scandir(args.sample_dir) if f.is_dir()]

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            sample_name = subfolder + calc.dir_sep + sample + "_AA_clean.txt"

            if os.path.isfile(sample_name) and os.path.getsize(sample_name) > 0:
                logging.debug("processing {}".format(sample_name))
                calc.run_calc(sample)

if __name__ == "__main__":
    run_calc(sys.argv[1:])

#TODO for testing, do not use in production
# sample_dir = r"D:\17 Dutihl Lab\_tools\diversitools"
# #run all samples
# run_calc(["-d", sample_dir, "-a"])

#or run one sample, or a list of
# sample = "MGXDB009139"
# run_calc(["-d", sample_dir, "-s", sample])
