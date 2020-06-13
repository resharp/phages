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
                              ,   usecols=[0, 1, 2, 3, 10]
                              ,   names=["Protein", "gene_fam", "gene_fam_annot", "region", "Annotation"]
                              )
        anno_df["ref"] = ref

        return anno_df

    def read_files(self):

        logging.debug("start reading tables")

        self.gene_df = self.read_gene_annotation("crassphage_refseq")

        logging.debug("end reading tables")

    def merge_files(self):

        data = self.gene_df

        merge_df = data
        # add annotation, e.g. region
        # merge_df = data.merge(self.gene_anno_df
        #                       , left_on=data.protein
        #                       , right_on=self.gene_anno_df.Protein
        #                       , how='left').drop(["key_0", "Protein", "ref_y"], axis=1)
        # merge_df.rename(columns={"ref_x": "ref"}, inplace=True)

        self.gene_df = merge_df

    def write_summary(self):

        filename = "{}{}out_gene_all_statistics.txt".format(
            self.sample_dir, self.dir_sep)
        self.gene_df.to_csv(path_or_buf=filename, sep='\t', index=False)

    def do_analysis(self):

        self.read_files()

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

