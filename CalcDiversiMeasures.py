import argparse
import logging
import os
import sys
import numpy as np
import pandas as pd

# we will calculate
# per genome AND per gene
#   average AA coverage
#   standard dev
#   normalized stdev
#   sum(CntNonSys)
#   sum(CntSyn)
#   dN/dS
#   pN/pS? -> search definition
#   % TopCodoncnt
#   % SndCodoncnt
#   % TrdCodoncnt
#   Length of gene
#
#   then filter based on threshold, e.g.
#       > 75% gene average (or a cut-off on normalized stddev)
#       > 75% genome average
#   first we do this on one sample and output a filtered file with measures per gene

class CalcDiversiMeasures:

    aa_table_name = ""
    codon_table_name = ""
    gene_table_name = ""

    sample_dir = ""
    sample = ""

    dir_sep = ""
    aa_df = None
    gene_df = None
    codon_df = None

    def __init__(self, sample_dir, sample):

        logging.basicConfig(filename='CalcDiversiMeasures.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.DEBUG)

        self.sample_dir = sample_dir
        self.sample = sample

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

    def read_files(self):
        logging.debug("start reading tables")

        self.read_aa_table()

        self.read_codon_table()

        self.integrate_tables()

        logging.debug("finished reading tables")

    def read_aa_table(self):

        self.aa_table_name = self.sample_dir + self.dir_sep + self.sample + self.dir_sep + self.sample + "_AA_clean.txt"
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

        #before aggregating first make a filtered derived measure for SndAAcnt_perc
        #we only want to keep percentages > 1% in the sample
        self.aa_df["SndAAcnt_perc"] = self.aa_df["SndAAcnt"] / self.aa_df["AAcoverage"]
        self.aa_df["CntSnp"] = self.aa_df["CntSyn"] + self.aa_df["CntNonSyn"]

        #TODO: filter out nan

        # debug_data[debug_data['SndAAcnt_perc'] > 0.01]

        self.aa_df['SndAAcnt_perc_filtered'] = 0
        self.aa_df.loc[self.aa_df["SndAAcnt_perc"] > 0.01, 'SndAAcnt_perc_filtered'] = self.aa_df["SndAAcnt_perc"]

        debug_data = self.aa_df[["Protein", "AAPosition", "SndAAcnt_perc", "SndAAcnt_perc_filtered"]]

        #mean coverage of all genes
        genome_coverage_mean = self.aa_df.AAcoverage.mean().round(decimals=2)

        #now first determine the gene table
        #then calculate average coverage and stdev
        #and output to new file

        # proteins = self.aa_df.Protein.unique()
        # for protein in proteins:
        #     print(protein)

        self.gene_df = self.aa_df.groupby("Protein").agg(
            {
                'AAcoverage': ["mean", "std", "sum"],
                'TopAAcnt': ["mean", "std", "sum"],
                'SndAAcnt': ["mean", "std", "sum"],
                'TrdAAcnt': ["mean", "std", "sum"],
                'CntNonSyn': 'sum',
                'CntSyn': 'sum',
                'SndAAcnt_perc': 'mean',
                'SndAAcnt_perc_filtered': 'mean',
                'syn': 'sum',
                'non_syn': 'sum',
                'CntSnp': 'median'
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

        #self.gene_df["SndAAcnt_perc"] = self.gene_df["SndAAcnt_sum"]/self.gene_df["AAcoverage_sum"]
        self.gene_df["dN/dS"] = self.gene_df["CntNonSyn_sum"]/self.gene_df["CntSyn_sum"]

        self.gene_df["syn_ratio"] = self.gene_df["syn_sum"] / self.gene_df["non_syn_sum"]

        self.gene_df["pN/pS"] = self.gene_df["syn_ratio"] * self.gene_df["CntNonSyn_sum"] / self.gene_df["CntSyn_sum"]

        np.seterr(divide='ignore')
        self.gene_df["log10_dN/dS"] = np.where(self.gene_df["dN/dS"] > 0, np.log10(self.gene_df["dN/dS"]), 0)
        self.gene_df["log10_pN/pS"] = np.where(self.gene_df["pN/pS"] > 0, np.log10(self.gene_df["pN/pS"]), 0)
        np.seterr(divide='warn')

        #quality measures
        self.gene_df["AAcoverage_perc"] = self.gene_df.AAcoverage_mean / genome_coverage_mean
        #coefficient of variation for a gene
        self.gene_df["AAcoverage_cv"] = self.gene_df.AAcoverage_std / self.gene_df.AAcoverage_mean

        #round all floats to four decimals
        self.gene_df = self.gene_df.round(decimals=4)

        #gene output table should contain sample (for later integrating over multiple samples)
        self.gene_df["sample"] = self.sample
        # print(self.gene_df.dtypes)

    def write_measures(self):

        self.gene_table_name = self.sample_dir + self.dir_sep + self.sample + self.dir_sep + self.sample + "_gene_measures.txt"
        #the gene name is in the index
        self.gene_df.to_csv(path_or_buf=self.gene_table_name, sep='\t')



def run_calc(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    parser.add_argument("-s", "--sample", dest="sample",
                        help="sample with Diversitools input", metavar="[sample}", required=True)

    args = parser.parse_args(args_in)

    print("Start running CalcDiversiMeasures")
    print("sample_dir: " + args.sample_dir)
    print("sample: " + args.sample)

    calc = CalcDiversiMeasures(args.sample_dir, args.sample)

    calc.read_files()

    calc.calc_measures()

    calc.write_measures()

if __name__ == "__main__":
    run_calc(sys.argv[1:])

#TODO for testing, do not use in production
# samples = ["MGXDB000864","MGXDB008660", "MGXDB009139", "MGXDB023930", "MGXDB022262"]
# #samples = ["MGXDB009139"]
#
# sample_dir = r"D:\17 Dutihl Lab\_tools\diversitools"
# for sample in samples:
#     run_calc(["-d", sample_dir, "-s", sample])
