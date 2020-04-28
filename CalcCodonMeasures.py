import argparse
import logging
import os
import pandas as pd
from scipy.stats import entropy
import sys

# CalcCodonMeasures aggregates on codon level over all samples based on the DiversiTools _AA.txt file
#
# output will be put in CodonMeasures (later to be processed by MakeDiversiPlots)
#
# to do: add region in the output
#
# We calculate per codon (by pile up/stacking of all reads of all samples on top of each other):
# - the overall coverage,
# - the total number of possible codon types for that position (may be 1)
# - and the entropy at that position (over all samples), [not including the ref genome codon type!]
# - added #syn and #non-syn to output
# - to do:
# (we can later aggregate this per gene, to have less data points).
# Goal is to have a trend line for the dependence of the entropy or SNP density (micro diversity)
#
# results will be put in CodonMeasures directory (sub dir of sample_dir)
# plots will be made by different class

# Also, we want to split the samples based on age_category, and therefore we need the sample metadata.
# First calculate the overall entropy, and later make a stack for the age categories.


class CalcCodonMeasures:

    sample_dir = ""
    ref = ""
    codon_table = ""

    aa_df = None
    sample_meta_df = None
    sample_measures_df = None
    codon_df = None

    out_dir = None

    def __init__(self, sample_dir, ref, codon_table):

        self.sample_dir = sample_dir
        self.ref = ref
        self.codon_table = codon_table

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'CalcCodonMeasures.log', filemode='w',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level=logging.DEBUG)

    def read_and_concat_measures(self, file_name_ext, ref=None, usecols = []):

        subfolders = [f.path for f in os.scandir(self.sample_dir)
                      if f.is_dir() and ref in f.name and "GenePlots" not in f.name]

        samples = []

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            sample = sample.split("_")[0]
            sample_name = subfolder + self.dir_sep + sample + file_name_ext

            # check if the sample is above breadth threshold based on sample measures
            filter_data = self.sample_measures_df[
                (self.sample_measures_df["sample"] == sample) & (self.sample_measures_df.ref == ref)].reset_index()
            if len(filter_data) == 0:
                logging.warning("No sample measures for {sample} and {ref}.".format(sample=sample, ref=ref))
                continue
            if filter_data.at[0, "breadth_1x"] < 0.05:
                logging.info("skipped sample {sample} for {ref}.".format(sample=sample, ref=ref))
                continue

            if os.path.isfile(sample_name) and os.path.getsize(sample_name) > 0:
                logging.info("processing {}".format(sample_name))

                if len(usecols) != 0:
                    measure_df = pd.read_csv(sample_name
                                             ,  sep='\t'
                                             , usecols=usecols
                                             )
                else:
                    measure_df = pd.read_csv(sample_name
                                             ,  sep='\t'
                                             )

                measure_df["sample"] = sample
                samples.append(measure_df)

        return pd.concat(samples)

    @staticmethod
    def age_category(row):
        return row.sample_name.split("_")[1]

    def read_sample_metadata(self):

        metadata_filename = self.sample_dir + self.dir_sep + "metadata_ERP005989.txt"

        self.sample_meta_df = pd.read_csv(metadata_filename
                                          , sep='\t'
                                          )

        self.sample_meta_df = self.sample_meta_df[['analysis.run.accession', 'sample.sample_name']]

        self.sample_meta_df.rename(columns={'analysis.run.accession': 'run',
                                            'sample.sample_name': 'sample_name'}, inplace=True)

        self.sample_meta_df["age_cat_short"] = self.sample_meta_df.apply(self.age_category, axis=1)

    def read_sample_measures(self):

        # sample measures contains the breadth_1x, ..10x, ..50x and ..100x fractions
        measure_file_name = self.sample_dir + self.dir_sep + "sample_measures.txt"
        self.sample_measures_df = pd.read_csv(measure_file_name,
                                              sep="\t")

    def read_codon_table(self):

        self.codon_df = pd.read_csv(self.codon_table
                                    ,   sep=","
                                    ,   usecols=[0,2,3,4])

    def read_files(self):

        logging.debug("start reading tables")

        self.read_sample_measures()

        self.read_sample_metadata()

        self.aa_df = self.read_and_concat_measures("_AA_clean.txt", self.ref,
                                                   [2, 3, 6, 10, 11, 19, 20, 21, 22, 23, 24, 25])

        self.read_codon_table()

        logging.debug("end reading tables")

    def merge_files(self):
        # add age_cat_short to self.aa_df
        self.aa_df = self.aa_df.merge(self.sample_meta_df[["run", "age_cat_short"]],
                                      left_on=self.aa_df["sample"],
                                      right_on=self.sample_meta_df["run"],
                                      how="inner")\
            .drop(["key_0", "run"], axis=1)
        self.aa_df.rename(columns={'age_cat_short': 'age_cat'}, inplace=True)

        self.aa_df = self.aa_df.merge(self.codon_df
                                    , left_on=self.aa_df.RefCodon
                                     , right_on=self.codon_df.codon
                                     , how='left')#.reset_index()

    def calc_and_write_measures(self):

        self.aa_df = self.aa_df[self.aa_df.TopCodon.notnull()]

        self.calc_and_write_measures_for_subset(self.aa_df, "all")

        age_cats = self.aa_df.age_cat.unique()
        for age_cat in age_cats:
            aa_df_age = self.aa_df[self.aa_df.age_cat == age_cat]

            # to do: write to one file instead of different files per reference genome
            self.calc_and_write_measures_for_subset(aa_df_age, age_cat)

    def calc_and_write_measures_for_subset(self, aa_df, age_cat=""):

        group = aa_df.groupby(["Protein", "AAPosition"])

        logging.debug("start pileup and calculation of entropy for age category {age_cat}".format(age_cat=age_cat))
        pileup_df = group.apply(self.count_codons)

        pileup_df["snp"] = 0
        pileup_df.loc[pileup_df.nr_snps > 0, "snp"] = 1

        pileup_df["age_cat"] = age_cat

        logging.debug("end pileup and calculation of entropy")

        filename = self.out_dir + self.dir_sep + "codon_entropy_{ref}_{age_cat}.txt"\
            .format(ref=self.ref, age_cat=age_cat)

        pileup_df.round(decimals=4).to_csv(path_or_buf=filename, sep='\t', index=False)

    def count_codons(self, group):

        # we want to calculate the entropy in a set-based manner
        # so first see what codons are unique
        # and then sum

        group_1 = pd.DataFrame({"codon": group["TopCodon"], "count": group["TopCodoncnt"]})
        group_2 = pd.DataFrame({"codon": group["SndCodon"], "count": group["SndCodoncnt"]})
        group_3 = pd.DataFrame({"codon": group["TrdCodon"], "count": group["TrdCodoncnt"]})

        codons = pd.concat([group_1, group_2, group_3])

        codon_sums = codons.groupby("codon").sum()

        aa_coverage_sum = group.AAcoverage.sum()

        # how many codons next to the most abundant have an abundance of at least 3 an more than 1% in population?
        possible_snps_including_major = codon_sums[(codon_sums["count"] > 2) &
                                                   (codon_sums["count"]/aa_coverage_sum > 0.01)]

        nr_snps = possible_snps_including_major.count()[0] - 1

        codon_entropy = entropy(codon_sums, base=10)

        df_return = pd.DataFrame({
            "protein": group["Protein"].head(1),
            "position": group["AAPosition"].head(1),
            "CntSyn": group["CntSyn"].head(1),
            "CntNonSyn": group["CntNonSyn"].head(1),
            "syn": group["syn"].head(1),
            "non_syn": group["non_syn"].head(1),
            "entropy": codon_entropy,
            "coverage": aa_coverage_sum,
            "nr_snps": nr_snps }
        )

        return df_return

    def create_output_dir(self):

        self.out_dir = self.sample_dir + self.dir_sep + "CodonMeasures"

        os.makedirs(self.out_dir, exist_ok=True)

    def run_calc(self):

        self.read_files()

        self.merge_files()

        self.create_output_dir()

        self.calc_and_write_measures()


def run_calc(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)

    parser.add_argument("-r", "--ref", dest="ref",
                        help="reference genome id", metavar="[ref]", required=True)

    parser.add_argument("-c", "--codon_table", dest="codon_table",
                        help="file with probabilities of syn or non syn mutations",
                        metavar="[codon_table]", required=True)

    args = parser.parse_args(args_in)

    print("Start running CalcCodonMeasures")
    print("sample_dir: " + args.sample_dir)
    print("reference genome: " + args.ref)

    calc = CalcCodonMeasures(args.sample_dir, args.ref, args.codon_table)

    calc.run_calc()


if __name__ == "__main__":
    run_calc(sys.argv[1:])

# TODO for testing, do not use in production
# sample_dir = r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"
# codon_table = r"D:\17 Dutihl Lab\source\phages_pycharm\input_files\codon_syn_non_syn_probabilities.txt"
# #
# # or run one sample
# ref = "crassphage_refseq"
# run_calc(["-d", sample_dir, "-r", ref, "-c", codon_table])