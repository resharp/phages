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
# we do this on one sample and output a file with measures per gene
# there is no filtering here
#
# the parameter ref (e.g. cs_ms_21) is only used for determination of read/write directory:
# [sample]_[ref] e.g. ERR526006_cs_ms_21
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

        self.sample_dir = sample_dir

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'CalcDiversiMeasures.log', filemode='w',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level=logging.DEBUG)

    def read_files(self, sample, ref):
        logging.debug("start reading tables")

        self.read_aa_table(sample, ref)

        self.read_codon_table()

        self.integrate_tables(sample)

        logging.debug("finished reading tables")

    def read_aa_table(self, sample, ref):

        self.aa_table_name = self.sample_ref_dir(sample, ref) + sample + "_AA_clean.txt"

        self.aa_df = pd.read_csv(self.aa_table_name
                                 ,  sep='\t'
                                 ,  usecols=range(2,26)
                                 )
        # automatic columns with first two skipped
        # SKIPPED: Sample	Chr
        # Protein	AAPosition	RefAA	RefSite	RefCodon	FstCodonPos	SndCodonPos	TrdCodonPos
        # CntNonSyn	CntSyn	NbStop	TopAA	TopAAcnt	SndAA	SndAAcnt	TrdAA	TrdAAcnt
        # TopCodon	TopCodoncnt	SndCodon	SndCodoncnt	TrdCodon	TrdCodoncnt	AAcoverage

    def sample_ref_dir(self, sample, ref):
        return self.sample_dir + self.dir_sep + sample + "_" + ref + self.dir_sep

    def read_codon_table(self):
        self.codon_table_name = self.sample_dir + self.dir_sep + "codon_syn_non_syn_probabilities.txt"

        self.codon_df = pd.read_csv(self.codon_table_name
                                    ,   sep=","
                                    ,   usecols=[0,2,3,4])

    def integrate_tables(self, sample):

        # sometimes we have mapped reads and they do not map in the coding regions
        # therefore, the column RefCodon is empty (as are all other columns) and we cannot join and should
        # stop processing to prevent further errors with merging of columns
        nr_unique_refcodons = self.aa_df.RefCodon.nunique()
        if np.isnan(nr_unique_refcodons) or nr_unique_refcodons < 2:
            error_message = "processing will be terminated to prevent errors. " + \
                            "not enough unique refcodons (={nr}) ".format(nr=nr_unique_refcodons) + \
                           "for sample {sample}".format(sample=sample)

            logging.error(error_message)
            raise Exception(error_message)

        else:
            info_message = "nr of unique refcodons {nr} ".format(nr=nr_unique_refcodons) + \
                           "for sample {sample}".format(sample=sample)
            logging.info(info_message)

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

        self.calculate_pn_ps_sliding_window(window_size=60)

    def calculate_pn_ps_sliding_window(self, window_size):

        # pN/pS calculation for sliding window
        # add simple moving averages for syn, non_syn, CntSyn and CntNonSyn
        ext = "_SMA_{}".format(window_size)

        self.aa_df["CntSyn" + ext] = self.aa_df["CntSyn"].replace(np.nan,0).rolling(window=window_size).mean()
        self.aa_df["CntNonSyn" + ext] = self.aa_df["CntNonSyn"].replace(np.nan,0).rolling(window=window_size).mean()
        self.aa_df["syn" + ext] = self.aa_df["syn"].replace(np.nan,0).rolling(window=window_size).mean()
        self.aa_df["non_syn" + ext] = self.aa_df["non_syn"].replace(np.nan,0).rolling(window=window_size).mean()

        # disregard values where AAPosition < window_size
        # otherwise you get a moving average partly based on values of the previous gene
        self.aa_df.loc[self.aa_df.AAPosition < window_size, "CntSyn" + ext] = 0
        self.aa_df.loc[self.aa_df.AAPosition < window_size, "CntNonSyn" + ext] = 0
        self.aa_df.loc[self.aa_df.AAPosition < window_size, "syn" + ext] = 0
        self.aa_df.loc[self.aa_df.AAPosition < window_size, "non_syn" + ext] = 0

        # median for whole sample
        count_snp_median = self.aa_df.CntSnp[self.aa_df["CntSnp"] != 0].median().round(decimals=4)
        snp_pseudo_count = np.sqrt(count_snp_median) / 2

        # this might result in infinity which is ok for now
        self.aa_df["syn_ratio" + ext] = \
            self.aa_df["syn" + ext].div(self.aa_df["non_syn" + ext])

        self.aa_df["pN_pS" + ext] = self.aa_df["syn_ratio" + ext] * \
            (self.aa_df["CntNonSyn" + ext] + snp_pseudo_count)/(self.aa_df["CntSyn" + ext] + snp_pseudo_count)

        self.aa_df["log10_pN_pS" + ext] = np.log10(self.aa_df["pN_pS" + ext])

        self.aa_df = self.aa_df.join(self.aa_df.groupby('Protein')['AAPosition'].max(), on='Protein', rsuffix='_max')

        shift_back = int(window_size / 2)

        measure = "log10_pN_pS_{}".format(window_size)
        self.aa_df[measure] = \
            self.aa_df["log10_pN_pS" + ext].shift(periods=-shift_back)

        # after shift make right side empty in order to disregard values of next gene
        self.aa_df.loc[self.aa_df.AAPosition_max - self.aa_df.AAPosition < shift_back, measure] = np.nan

        # round all values to four decimals
        self.aa_df = self.aa_df.round(decimals=4)

    def create_plot_dir(self, sample, ref):

        self.plot_dir = self.sample_ref_dir(sample, ref) + "DiversiMeasures"
        os.makedirs(self.plot_dir, exist_ok=True)

    # aggregate measures
    def aggregate_measures(self, sample):

        # mean coverage of all genes

        # TODO add CntSnp median over all genes
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

        # remove multi-index set on column axis
        # https://www.shanelynn.ie/summarising-aggregation-and-grouping-data-in-python-pandas/
        self.gene_df.columns = ["_".join(x) for x in self.gene_df.columns.ravel()]

        self.gene_df.AAcoverage_sum = pd.to_numeric(self.gene_df.AAcoverage_sum, downcast='unsigned', errors='coerce')
        self.gene_df.CntNonSyn_sum = pd.to_numeric(self.gene_df.CntNonSyn_sum, downcast='unsigned', errors='coerce')
        self.gene_df.CntSyn_sum = pd.to_numeric(self.gene_df.CntSyn_sum, downcast='unsigned', errors='coerce')
        self.gene_df.TopAAcnt_sum = pd.to_numeric(self.gene_df.TopAAcnt_sum, downcast='unsigned', errors='coerce')
        self.gene_df.SndAAcnt_sum = pd.to_numeric(self.gene_df.SndAAcnt_sum, downcast='unsigned', errors='coerce')
        self.gene_df.TrdAAcnt_sum = pd.to_numeric(self.gene_df.TrdAAcnt_sum, downcast='unsigned', errors='coerce')

        # derived measures
        self.gene_df["syn_ratio"] = self.gene_df["syn_sum"] / self.gene_df["non_syn_sum"]

        # take median of non-zero values over all amino acids (alle genes)
        count_snp_median = self.aa_df.CntSnp[self.aa_df["CntSnp"] != 0].median().round(decimals=4)
        snp_pseudo_count = np.sqrt(count_snp_median) / 2

        self.gene_df["pN/pS"] = self.gene_df["syn_ratio"] * \
            (self.gene_df["CntNonSyn_sum"] + snp_pseudo_count)/(self.gene_df["CntSyn_sum"] + snp_pseudo_count)

        self.gene_df["log10_pN/pS"] = np.where(self.gene_df["pN/pS"] > 0, np.log10(self.gene_df["pN/pS"]), 0)

        # quality measures
        self.gene_df["AAcoverage_perc"] = self.gene_df.AAcoverage_mean / genome_coverage_mean
        # coefficient of variation for a gene
        self.gene_df["AAcoverage_cv"] = self.gene_df.AAcoverage_std / self.gene_df.AAcoverage_mean

        # round all floats to four decimals
        self.gene_df = self.gene_df.round(decimals=4)

        # gene output table should contain sample (for later integrating over multiple samples)
        self.gene_df["sample"] = sample

    def write_aggregated_measures(self, sample, ref):

        self.gene_table_name = self.sample_ref_dir(sample, ref) + sample + "_gene_measures.txt"
        # the gene name is in the index
        self.gene_df.to_csv(path_or_buf=self.gene_table_name, sep='\t')

    def write_bin_measures(self, sample, ref):

        bin_measure_name = self.sample_ref_dir(sample, ref) + sample + "_bin_measures.txt"

        self.aa_df['sample'] = sample
        data = self.aa_df[['sample', 'Protein', 'AAPosition', 'AAcoverage', 'entropy', 'log10_pN_pS_60']]

        # TODO: also join self.gene_df and add
        # AAcoverage_perc	AAcoverage_cv
        # in order to filter on it in further data integration steps
        data = data.join(self.gene_df[["AAcoverage_perc", "AAcoverage_cv"]],
                  on='Protein', rsuffix='_gene')

        data.to_csv(path_or_buf=bin_measure_name, sep='\t', index=False)

    # now here is a region with verbose code that is not part of the main data analysis
    #region verbose
    def line_plots_for_coverage_and_pn_ps(self):

        self.line_plots_for_measure("AAcoverage", ylim_bottom=0)
        self.line_plots_for_measure("log10_pN_pS_60")

    def line_plots_for_measure(self, measure, ylim_bottom=None):

        data = self.aa_df[["Protein", "AAPosition", "AAcoverage", "log10_pN_pS_60", "AAPosition_max"]]
        genes = data.Protein.unique()

        # TODO: determine the number of genes and make a plot for every set of 6 genes
        nr_plots = int(np.round((len(genes)/6)))

        # we loop through the plots and subplots and then select the next gene instead of looping through the genes
        i = 0
        # we always make one extra plot (so some of the sub plots of the last plot may be empty)
        for j in range(0,nr_plots):

            fig, axs = plt.subplots(3,2)
            fig.tight_layout(pad=2)
            for row in axs:
                for ax in row:
                    # to prevent
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

    def joint_plot_for_aa_changes_and_entropy(self):
        # here we show the AAPosition the x-axis and the entropy on the y-axis for every SCP
        # though entropy seems to be an attribute of a sample we would like to see if there is

        self.aa_df["SCP"] = 0
        self.aa_df.loc[self.aa_df.SndAAcnt_perc_polymorphism > 0.01, "SCP"] = 1

        data = self.aa_df[["Protein", "AAPosition", "TopCodon", "RefCodon", "SndAAcnt_perc_polymorphism", "SCP", "entropy"]]

        # to do: do this for all genes in verbose mode?
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

                figure_name = self.plot_dir + self.dir_sep + "joint_plot_{gene}.png".format(gene=gene)

                plt.savefig(figure_name)

    def find_and_write_local_regions(self, sample, ref):

        # now we want to determine the distances between SNPs in self.aa_df
        # and then take the 5 percentile shortest distances to identify regions that are close
        # for now we determine SNP as non-empty RefCodon unequal to non-empty TopCodon

        # SCP "single codon polymorphism" = difference between TopCodon and RefCodon
        scp_df = self.aa_df
        scp_df = scp_df[scp_df.TopCodon.notnull()]
        scp_df = scp_df[scp_df.RefCodon.notnull()]
        scp_df = scp_df[scp_df.TopCodon != scp_df.RefCodon]

        self.find_and_write_hypervariable_regions(sample, ref, scp_df, "non syn regions compared to ref")

        # TODO: reconsider, finding high peaks of entropy might be interesting
        # however, it will not find broad ranges of high entropy (we have to find some balance for the height_percentile
        # and the perceptile in find_and_write_hypervariable_regions depending on our question)
        height_percentile = 95
        entropy_cutoff = np.percentile(self.aa_df[self.aa_df["entropy"].notnull()].entropy, height_percentile)

        entropy_df = self.aa_df[["Protein", "AAPosition", "SndAAcnt_perc_polymorphism", "entropy"]]
        entropy_df = entropy_df[(entropy_df.entropy > entropy_cutoff)]

        self.find_and_write_hypervariable_regions(sample, ref, entropy_df, "high entropy regions")

        # we should also be able to do the same for any other positions
        # for example for distances between within-sample variations

    # poi_df contains Protein and positions of interest of filtered rows
    # poi_df only only needs to contain Protein and AAPosition and you can put anything in AAPosition
    # however, for debugging purposes, you might want to add filter fields like RefCodon and TopCodon
    # that have contributed to the determination of the positions of interest
    def find_and_write_hypervariable_regions(self, sample, ref, poi_df, title):

        # poi_df = poi_df[["Protein", "AAPosition", "RefCodon", "TopCodon"]].reset_index()
        poi_df = poi_df[["Protein", "AAPosition"]].reset_index()

        #add some extra fields to use for positioning
        poi_df.rename(columns={'index': 'AAPositionGenome'}, inplace=True)
        poi_df['AAPositionGenome'] = poi_df['AAPositionGenome'] + 1 #+1 to correct for 0-based result of reset_index()

        poi_df = poi_df.reset_index()
        poi_df.rename(columns={'index': 'PoiPosition'}, inplace=True)
        poi_df['PoiPosition'] = poi_df['PoiPosition'] + 1 #+1 to correct for 0-based result of reset_index()

        # distance_bw is the backward distance
        poi_df['distance_bw'] = poi_df['AAPositionGenome'].diff()
        # distance_fw is the forward distance
        poi_df['distance_fw'] = poi_df['distance_bw'].shift(-1)

        # these next two lines are the core of the method: determine the 5% percentile
        # of the distance distribution between positions of interest
        percentile = 5
        distance_cutoff = np.percentile(poi_df[poi_df['distance_bw'].notnull()].distance_bw, percentile)

        # for experimenting with manually set distance cut-offs, not determined by percentile
        # TODO ************************
        # distance_cutoff = 1

        # apply filter on distance cutoff between POIs (positions-of-interest)
        poi_df = poi_df[(poi_df.distance_bw <= distance_cutoff) | (poi_df.distance_fw <= distance_cutoff)]

        # PoiDistance = the distance between FILTERED rows
        # =>PoiDistance 1 means that these SCPs should be joined in a hypervariable region according to distance cut-off
        # e.g. when the cut-off is set at 2 the "distance" between rows can be 2, and still the PoiDistance is 1
        # first calculate the backward distance
        poi_df['PoiDistance_bw'] = poi_df['PoiPosition'].diff()
        # forward distance is just the backward distance of the next row
        poi_df['PoiDistance_fw'] = poi_df['PoiDistance_bw'].shift(-1)

        # now filter out just the joining regions
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
            # initialize previous position in such a way that the first condition will hold and a region will be started
            old_pos = positions[0] - distance_cutoff
            # always start a new region for each gene
            region = []
            for pos in positions:
                if (pos <= old_pos + distance_cutoff):
                    region.append(pos)
                else:
                    # first finish the region
                    if len(region) > 0:
                        first_aa_pos = region[0]
                        last_aa_pos = region[-1]
                        data_for_df.append([row["Protein"], first_aa_pos, last_aa_pos])
                    # then start a new region
                    region = []
                    region.append(pos)
                # save previous position
                old_pos = pos
            # finish the last region
            if len(region) > 0:
                first_aa_pos = region[0]
                last_aa_pos = region[-1]
                data_for_df.append([row["Protein"], first_aa_pos, last_aa_pos])

        columns = ["Protein", "AAPositionStart", "AAPositionEnd"]

        new_gene_poi_df = pd.DataFrame(columns=columns, data=data_for_df)

        new_gene_poi_df["length"] = new_gene_poi_df["AAPositionEnd"] - new_gene_poi_df["AAPositionStart"] + 1

        gene_region_table_name = self.sample_ref_dir(sample, ref) + \
                                 sample + "_gene_" + title.replace(" ", "_") + ".txt"
        new_gene_poi_df.to_csv(gene_region_table_name, sep='\t', index=False)

    #endregion verbose

    def run_calc(self, sample, ref):

        self.read_files(sample, ref)

        self.calc_measures()

        self.create_plot_dir(sample, ref)

        verbose = False
        if verbose:
            self.joint_plot_for_aa_changes_and_entropy()
            self.line_plots_for_coverage_and_pn_ps()

            # this is still too buggy for the current pipeline (when there are not so many reads mapped)
            self.find_and_write_local_regions(sample, ref)

        self.aggregate_measures(sample)

        self.write_aggregated_measures(sample, ref)

        self.write_bin_measures(sample, ref)


def run_calc(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    parser.add_argument("-s", "--sample", dest="sample",
                        help="sample with Diversitools input", metavar="[sample}", required=False)

    parser.add_argument("-r", "--ref", dest="ref",
                        help="reference genome id", metavar="[ref}", required=True)

    parser.add_argument("-a", "--all", action="store_true", default=False,
                        help="run all samples in sample directory", required=False)


    args = parser.parse_args(args_in)

    print("Start running CalcDiversiMeasures")
    print("sample_dir: " + args.sample_dir)
    print("reference genome: " + args.ref)
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
        calc.run_calc(args.sample, args.ref)
        return

    if args.all:
        subfolders = [f.path for f in os.scandir(args.sample_dir) if f.is_dir()
                      and args.ref in f.name and not ("GenePlots" in f.name)]

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            sample = sample.split("_")[0]
            sample_name = subfolder + calc.dir_sep + sample + "_AA_clean.txt"

            if os.path.isfile(sample_name) and os.path.getsize(sample_name) > 0:
                logging.debug("processing {}".format(sample_name))
                calc.run_calc(sample, args.ref)

if __name__ == "__main__":
    run_calc(sys.argv[1:])

#TODO for testing, do not use in production
# sample_dir = r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"
# ref = "crassphage_refseq"
# #run all samples
# run_calc(["-d", sample_dir, "-r", ref, "-a"])

#or run one sample, or a list of
# sample = "ERR525804"
# run_calc(["-d", sample_dir, "-s", sample, "-r", ref])
