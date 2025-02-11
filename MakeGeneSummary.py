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
#   	crassphage_refseq_gene_list
# T1	gene_fam_sample.0.95.10x.txt
# T2	crassphage_pN_pS_values.0.95.10x.txt
# T3	across_sample_gene_measures_table.txt
# hmm_hits	hmm_hits for copy numbers (first add genus?)
# 	ids_ref_genomes.txt
# for this we need a complete list with protein -> genus (not just for which we have measures)
# in the family analysis we read all separate
# with read_all_annotations

class MakeGeneSummary:

    sample_dir = ""

    gene_df = None
    genus1_stats_df = None
    fam_sample_df = None
    fam_df = None
    hhm_hits_df = None
    anno_df = None
    ref_meta_df = None
    pangenome_matrix_df = None
    across_sample_measures_df = None
    fam_snp_density_df = None

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

    def read_across_sample_measures(self):

        stats_filename = self.sample_dir + self.dir_sep + "filtered_across_sample_gene_measures_table.txt"

        df = pd.read_csv(stats_filename
                              ,   sep='\t'
                              ,   usecols=[0,1,12,17]
                              )

        df = df[df.age_cat == "all"]

        return df

    def read_hmm_hits(self):

        filename = self.sample_dir + self.dir_sep + "hmm_hits.tsv"

        df = pd.read_csv(filename
                         ,  sep='\t'
                         ,  names=["gene", "gene_fam", "e_value"]
                         ,  usecols=[0,1,2]
                         )
        logging.debug("Nr of HMM hits: " + str(len(df)))

        df = df[df.e_value < 1e-10]

        logging.debug("Nr of HMM hits with e-value < 1e-10: " + str(len(df)))

        # to do: return to 1e-30 threshold?
        # df = df[df.e_value < 1e-30]
        # logging.debug("Nr of HMM hits with e-value < 1e-30: " + str(len(df)))

        return df

    def read_all_annotations(self):

        ref_dir = self.sample_dir + self.dir_sep + "gene_list"

        files = [f.path for f in os.scandir(ref_dir) if not f.is_dir() and "_gene_list.txt" in f.name]

        anno_list = []

        for file in files:

            anno_df = pd.read_csv(file
                                  , sep='\t'
                                  , header=None
                                  , usecols=[0]
                                  , names=["Protein"]
                                  , skiprows=1
                                  )
            ref = file.split(self.dir_sep)[-1].replace("_gene_list.txt", "")
            anno_df["ref"] = ref

            anno_list.append(anno_df)

        df = pd.concat(anno_list)

        self.read_ref_metadata()

        ref_file_name = self.sample_dir + self.dir_sep + "ref_genome_ids.txt"

        df_ref = pd.read_csv(ref_file_name
                             , sep="\t"
                             , header=None
                             , comment="#"
                             , names=["ref", "genus"])
        df_ref.genus = df_ref.apply(self.shorten_genus, axis=1)

        df = df.merge(df_ref,
                      left_on=df.ref,
                      right_on=df_ref.ref,
                      how="inner").drop(["key_0", "ref_x"], axis=1)

        df.rename(columns={'ref_y': 'ref'}, inplace=True)
        return df

    def read_ref_metadata(self):

        ref_file_name = self.sample_dir + self.dir_sep + "ref_genome_ids.txt"

        self.ref_meta_df = pd.read_csv(ref_file_name
                                       , sep="\t"
                                       , header=None
                                       , comment="#"
                                       , names=["ref", "genus"]
                                       )
        self.ref_meta_df.genus = self.ref_meta_df.apply(self.shorten_genus, axis=1)

    @staticmethod
    def shorten_genus(row):
        return row.genus.replace("genus_", "")

    def read_files(self):

        logging.debug("start reading tables")

        self.gene_df = self.read_gene_annotation("crassphage_refseq")

        self.genus1_stats_df = self.read_genus1_stats()

        self.fam_sample_df = self.read_fam_sample_stats()

        self.across_sample_measures_df = self.read_across_sample_measures()

        self.hhm_hits_df = self.read_hmm_hits()

        self.anno_df = self.read_all_annotations()

        logging.debug("end reading tables")

    def aggregate_gene_fam(self):

        data = self.fam_sample_df

        data = data[["gene_fam", "log10_pN/pS_mean", "log10_pN/pS_count"]]

        # results are already aggregated, drop duplicates and be left with a single line per gene family
        df = data.drop_duplicates()

        self.fam_df = df

        df_snp = self.across_sample_measures_df.groupby(["gene_fam"]).agg(
            {'snp_density': 'mean'}
        ).reset_index()

        df_snp.rename(columns={'snp_density': 'snp_density_fam'}, inplace=True)

        self.fam_snp_density_df = df_snp

    def do_pangenome_analysis(self):

        data = self.hhm_hits_df

        # first add genus
        data = data.merge(self.anno_df,
                          left_on=data.gene,
                          right_on=self.anno_df.Protein,
                          how="inner").drop(["key_0", "Protein"], axis=1)

        # aggregate on a combination of gene_fam en genus and count the number of proteins
        df = data.groupby(["gene_fam", "ref"]).agg(
            {'gene': 'count'}
        ).reset_index()

        df.rename(columns={'gene': 'genus_count'}, inplace=True)

        # to do: first aggregate on gene_fam and count nr genera
        count_df = df.groupby(["gene_fam"]).agg(
            {'genus_count': 'count'}
        ).reset_index()

        df = df.set_index(["gene_fam", "ref"])

        # manual correction: we know that genus 4 has 1 portal protein
        df.loc["portal", "srr4295175_ms_5"] = 1

        # convert from multi-index to matrix
        df = df.unstack()
        df.columns = ["_".join(x) for x in df.columns.ravel()]
        df.columns = [x.replace("genus_count_","") for x in df.columns]
        df = df.reset_index()

        df = count_df.merge(df,
                            left_on=count_df.gene_fam,
                            right_on=df.gene_fam,
                            how='inner').drop(["key_0","gene_fam_y"], axis=1)
        df.rename(columns={'gene_fam_x': 'gene_fam'}, inplace=True)

        return df

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

        # add SNP densities on both gene and family level
        df_protein_snp_density = self.across_sample_measures_df[["protein", "snp_density"]]
        merge_df = merge_df.merge(df_protein_snp_density
                                  , left_on=merge_df.protein
                                  , right_on=df_protein_snp_density.protein
                                  , how="left").drop(["key_0", "protein_y"], axis=1)
        merge_df.rename(columns={"protein_x": "protein",
                                 "snp_density": "snp_density_gene"}, inplace=True)

        merge_df = merge_df.merge(self.fam_snp_density_df
                                  , left_on=merge_df.gene_fam
                                  , right_on=self.fam_snp_density_df.gene_fam
                                  , how="left").drop(["key_0", "gene_fam_y"], axis=1)
        merge_df.rename(columns={"gene_fam_x": "gene_fam"}, inplace=True)

        merge_df = merge_df.merge(self.pangenome_matrix_df
                                  , left_on=merge_df.gene_fam
                                  , right_on=self.pangenome_matrix_df.gene_fam
                                  , how='left').drop(["key_0", "gene_fam_y"], axis=1)
        merge_df.rename(columns={"gene_fam_x": "gene_fam"}, inplace=True)

        self.gene_df = merge_df

    def write_summary(self):

        filename = "{}{}out_gene_all_statistics.txt".format(
            self.sample_dir, self.dir_sep)
        self.gene_df.to_csv(path_or_buf=filename, sep='\t', index=False)

    def do_analysis(self):

        self.read_files()

        self.aggregate_gene_fam()

        self.pangenome_matrix_df = self.do_pangenome_analysis()

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

