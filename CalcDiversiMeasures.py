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
    gene_table_name = ""

    sample_dir = ""
    sample = ""

    dir_sep = ""
    aa_df = None
    gene_df = None

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

    def calc_measures(self):

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
                'CntSyn': 'sum'
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
        # self.gene_df["SndAAcnt_perc"] = self.gene_df["SndAAcnt_sum"]/self.gene_df["AAcoverage_sum"]
        self.gene_df["dN/dS"] = self.gene_df["CntNonSyn_sum"]/self.gene_df["CntSyn_sum"]

        self.gene_df["log10_dN/dS"] = np.log10(self.gene_df["dN/dS"])

        self.gene_df["AAcoverage_perc"] = self.gene_df.AAcoverage_mean / genome_coverage_mean
        #coefficient of variation for a gene
        self.gene_df["AAcoverage_cv"] = self.gene_df.AAcoverage_std / self.gene_df.AAcoverage_mean

        #round all floats to two decimals
        self.gene_df = self.gene_df.round(decimals=2)

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
# sample = "MGXDB008660"
# sample_dir = r"D:\17 Dutihl Lab\_tools\diversitools"
# run_calc(["-d", sample_dir, "-s", sample])
