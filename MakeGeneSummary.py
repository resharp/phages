import argparse
import logging
import os
import sys
import pandas as pd
import seaborn as sns;sns.set()

# MakeGeneSummary
#
# collect all stats for genes and (optionally) accompanying gene families in one export file
#
# input files
# T0	crassphage_refseq_gene_list
# T1	gene_fam_sample.0.95.10x.txt
# T2	crassphage_pN_pS_values.0.95.10x.txt
# T3	across_sample_gene_measures_table.txt
# hmm_hits	hmm_hits for copy numbers (first add genus?)
# 	ids_ref_genomes.txt


class MakeGeneSummary:

    sample_dir = ""

    gene_df = None
    genus1_stats_df = None
    fam_sample_df = None
    fam_df = None

    def __init__(self, sample_dir):

        self.sample_dir = sample_dir

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'MakeGeneSummary.log', filemode='w',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level=logging.DEBUG)

    def read_gene_annotation(self, ref):

        gene_anno_file_name = self.sample_dir + self.dir_sep + "{ref}_gene_list.txt".format(ref=ref)

        anno_df = pd.read_csv(gene_anno_file_name
                              ,   sep='\t'
                              ,   header=None
                              ,   skiprows=[0]
                              ,   usecols=[0, 1, 3, 10]
                              ,   names=["Protein", "gene_fam", "region", "Annotation"]
                              )
        return anno_df

    def read_genus1_stats(self):

        stats_filename = self.sample_dir + self.dir_sep + "crassphage_pN_pS_values.0.95.10x.txt"

        df = pd.read_csv(stats_filename
                              ,   sep='\t'
                              ,   usecols=[0, 3, 4, 5]
                              )
        return df

    def read_fam_sample_stats(self):

        stats_filename = self.sample_dir + self.dir_sep + "gene_fam_sample.0.95.10x.txt"

        df = pd.read_csv(stats_filename
                              ,   sep='\t'
                              ,   usecols=[1, 27, 33, 36, 41, 42, 45, 46]
                              )
        return df

    def read_files(self):

        logging.debug("start reading tables")

        self.gene_df = self.read_gene_annotation("crassphage_refseq")

        self.genus1_stats_df = self.read_genus1_stats()

        self.fam_sample_df = self.read_fam_sample_stats()

        logging.debug("end reading tables")

    def aggregate_gene_fam(self):

        data = self.fam_sample_df

        data = data[["gene_fam", "log10_pN/pS_mean", "log10_pN/pS_count"]]

        # results are already aggregated, drop duplicates and be left with a single line per gene family
        df = data.drop_duplicates()

        self.fam_df = df

    def merge_files(self):

        data = self.gene_df

        merge_df = data.merge(self.genus1_stats_df
                              , left_on=data.Protein
                              , right_on=self.genus1_stats_df.Protein
                              , how='left').drop(["key_0", "Protein_y"], axis=1)
        merge_df.rename(columns={"Protein_x": "protein",
                                 "log10_pN/pS_mean": "log10_pN/pS_mean_gene",
                                 "log10_pN/pS_count": "gene_in_samples"}, inplace=True)
        merge_df = merge_df.merge(self.fam_df
                                  , left_on=merge_df.gene_fam
                                  , right_on=self.fam_df.gene_fam
                                  , how='left').drop(["key_0", "gene_fam_y"], axis=1)
        merge_df.rename(columns={"gene_fam_x": "gene_fam",
                                 "log10_pN/pS_mean": "log10_pN/pS_mean_fam",
                                 "log10_pN/pS_count": "fam_in_samples"}, inplace=True)

        self.gene_df = merge_df

    def write_summary(self):

        filename = "{}{}out_gene_all_statistics.txt".format(
            self.sample_dir, self.dir_sep)
        self.gene_df.to_csv(path_or_buf=filename, sep='\t', index=False)

    def do_analysis(self):

        self.read_files()

        self.aggregate_gene_fam()

        self.merge_files()

        self.write_summary()


def do_analysis(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    args = parser.parse_args(args_in)

    make = MakeGeneSummary(args.sample_dir)

    make.do_analysis()


# if __name__ == "__main__":
#     run_calc(sys.argv[1:])


# TODO for testing, do not use in production
sample_dir = r"D:\17 Dutihl Lab\_tools\_summary"

do_analysis(["-d", sample_dir])

