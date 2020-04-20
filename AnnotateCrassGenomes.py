import argparse
import os

import pandas as pd


# we need to read the files and join the information in the following order:
# from Guerin dir
# ref_genes.tsv (genome, gene)
# hmm_hits.tsv  (gene, gene_fam, hhm_eval, hmm_score, AA_length)
#
# from Yutin dir
# to do: does this annotation contain new information?
# yutin_conserved_genes.txt (annotation, MSA-file = gene_fam, yutin_nr)
# crassphage_gene_list_add_yutin.txt    (protein, annotation, yutin_nr, region)


class AnnotateCrassGenomes:

    annotation_dir = ""
    genome_dir = ""

    ref_genes_name = "ref_genes.tsv"
    hmm_hits_name = "hmm_hits.tsv"

    yutin_genes_name = "yutin_conserved_genes.txt"
    crass_genes_name = "crassphage_gene_list_add_yutin.txt"

    pvog_annotation_name = "pvogs_annotations.tsv"
    pvog_hmm_results_name = "all_refs_pvog_table.txt"

    dir_sep = ""

    ref_genes_df = None
    hmm_hits_df = None
    yutin_genes_df = None
    crass_genes_df = None
    merge_df = None

    pvog_annot_df = None
    pvog_hmm_df = None

    def __init__(self, genome_dir, annotation_dir):

        self.annotation_dir = annotation_dir
        self.genome_dir = genome_dir

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

    def annotate(self):

        print("genomes in: " + self.genome_dir)
        print("annotation in: " + self.annotation_dir)

        self.read_files()

        self.join_files()

        self.write_files()

    def read_files(self):
        self.ref_genes_name = self.genome_dir + self.dir_sep + self.ref_genes_name
        self.hmm_hits_name = self.genome_dir + self.dir_sep + self.hmm_hits_name

        self.yutin_genes_name = self.annotation_dir + self.dir_sep + self.yutin_genes_name
        self.crass_genes_name = self.annotation_dir + self.dir_sep + self.crass_genes_name

        print("starting to read files")
        print(self.hmm_hits_name)
        print(self.ref_genes_name)
        print(self.yutin_genes_name)
        print(self.crass_genes_name)

        self.ref_genes_df = pd.read_csv(self.ref_genes_name
                                , header=None
                                , sep="\t"
                                , names=["genome","gene","meta"])
        print("number of genes in ref genomes: {nr_genes}".format(nr_genes=len(self.ref_genes_df)))

        self.hmm_hits_df = pd.read_csv(self.hmm_hits_name
                                       , header=None
                                       , sep="\t"
                                       , names=["gene", "gene_fam", "hmm_evalue", "hmm_score"]
                                       , usecols=[0, 1, 2, 3]
                                       #          , "start","end","AA_length"]
                                       # , dtype={"start": "int32", "end": "int32", "AA_length": "int32"}
                                       )

        # print(self.hmm_hits_df.dtypes)
        print("number of hmm hits: {nr_hits}".format(nr_hits=len(self.hmm_hits_df)))

        self.yutin_genes_df = pd.read_csv(self.yutin_genes_name
                                          , sep="\t"
                                          , index_col=None
                                          , usecols=range(0,23)
                                          )
        self.yutin_genes_df.rename(columns={'gene_annot': 'annot_yutin_hmm'}, inplace=True)

        # print(self.yutin_genes_df.dtypes)
        print("number of conserved genes: {nr_genes}".format(nr_genes=len(self.yutin_genes_df)))

        self.crass_genes_df= pd.read_csv(self.crass_genes_name
                                         , header=None
                                         , sep='\t'
                                         , index_col=None
                                         , usecols=[0,1,2,4,6,7,8,9,10]
                                         , names=["protein","start","end","annotation"
                                                    ,"yutin_gene_nr","region","occurs_in_fam","my_annot","remark2"]
                                         , dtype={"yutin_gene_nr": "str"})
        # we had to change yutin_gene_nr to dtype object in order to join it with yutin_genes_df.yutin_gene_nr (e.g.46N)
        print("number of proteins in annotated crAssphage: {nr_proteins}".format(nr_proteins=len(self.crass_genes_df)))
        self.crass_genes_df.rename(columns={'annotation': 'annot_genbank'}, inplace=True)

        self.read_pvog_files()

    def read_pvog_files(self):

        self.pvog_annotation_name = self.genome_dir + self.dir_sep + self.pvog_annotation_name

        # read the (processed) annotation from Nikos
        self.pvog_annot_df= pd.read_csv(self.pvog_annotation_name
                                        , sep="\t"
                                        )
        self.pvog_annot_df.rename(columns={'annotation_processed': 'annot_pvog'
                                           , 'annotation_raw': 'annot_pvog_raw'}, inplace=True)

        self.pvog_hmm_results_name = self.genome_dir + self.dir_sep + self.pvog_hmm_results_name

        # read the pvog hmm results
        self.pvog_hmm_df = pd.read_csv(self.pvog_hmm_results_name
                                       , delim_whitespace=True
                                       , comment="#"
                                       , usecols=[0, 2, 4, 5]
                                       , header=None
                                       , names=["gene", "pvog", "e_value", "score"]
                                       )

        # to do: what cut-off should we use here? 1e-30 or higher?
        # what happens if we take a more relaxed cut-off?
        # more relaxed cut-off is better for comparison with Guerin
        # self.pvog_hmm_df = self.pvog_hmm_df[self.pvog_hmm_df.e_value < 1e-10]
        self.pvog_hmm_df = self.pvog_hmm_df[self.pvog_hmm_df.e_value < 0.2]

        merge_df = self.pvog_hmm_df.merge(self.pvog_annot_df,
                                          left_on=self.pvog_hmm_df.pvog,
                                          right_on=self.pvog_annot_df.pvog,
                                          how="inner").drop(["key_0", "pvog_y"], axis=1)

        # sometimes we have two hits in pvog_hmm_df, we should take the best one
        # e.g. hvcf_a6_ms_4_36:
        # only take the best pvog hmm hit by sorting and then dropping gene duplicates
        merge_df = merge_df.sort_values(["gene", "score"], ascending=False)
        self.pvog_hmm_df = merge_df.drop_duplicates('gene')

        # debug = self.pvog_hmm_df[self.pvog_hmm_df.gene.str.match("hvcf_a6_ms_4")]

    def join_files(self):

        print("---------------")

        nr_fams_in_hits = len(self.hmm_hits_df.gene_fam.unique())
        print("number of fams in hits: {nr_fams}".format(nr_fams=nr_fams_in_hits))

        filtered_hits = self.hmm_hits_df[self.hmm_hits_df.hmm_evalue < 1e-30]

        print("number of hmm hits smaller than 1e-30: {nr_hits}".format(nr_hits=len(filtered_hits)))

        filtered_fams = self.yutin_genes_df[self.yutin_genes_df.gene_fam != " -"]
        print("number of gene fams with hmm profile: {nr_fams}".format(nr_fams=len(filtered_fams)))

        # now we start with all genes
        merge_df = self.ref_genes_df.merge(filtered_hits,
                                left_on=self.ref_genes_df.gene,
                                right_on=filtered_hits.gene,
                                how="left").drop(["key_0","meta", "gene_y"], axis=1)
        merge_df.rename(columns={'gene_x': 'gene'}, inplace=True)

        print("check nr of lines in merge_df: {nr_lines}".format(nr_lines=len(merge_df)))
        # merge_df.rename(columns={'gene_x': 'gene'}, inplace=True)

        merge_df = merge_df.merge(self.yutin_genes_df[["gene_fam","yutin_gene_nr","annot_yutin_hmm"]],
                                  left_on=merge_df.gene_fam,
                                  right_on=self.yutin_genes_df.gene_fam,
                                  how="left").drop(["key_0"], axis=1)

        # to do: consider joining on gene_fam instead of yutin gene number
        merge_df = merge_df.merge(self.crass_genes_df,
                                  left_on=merge_df.yutin_gene_nr,
                                  right_on=self.crass_genes_df.yutin_gene_nr,
                                  how="left").drop(["key_0"], axis=1)
        merge_df.rename(columns={'gene_fam_x': 'gene_fam'}, inplace=True)

        self.merge_df = merge_df

    def write_files(self):

        out_table_name = self.genome_dir + self.dir_sep + "out_gene_annotations.txt"

        self.merge_df.to_csv(path_or_buf=out_table_name, sep='\t', index=False)

        # now also write one file for every ref genome with only the annotations
        # with no annotation, make it "unknown"
        genomes = self.merge_df.genome.unique()

        list_dfs = []

        for genome in genomes:

            gene_list_name = self.genome_dir + self.dir_sep + "{genome}_gene_list.txt".format(genome=genome)
            genes_df = self.merge_df[
                self.merge_df.genome == genome][["gene", "gene_fam", "annot_yutin_hmm", "region", "annot_genbank"]]

            # if genome == 'crassphage_refseq', take region from original file (and not by joining on gene_fam)
            # it is already in self.crass_genes_df!
            if genome == "crassphage_refseq":
                genes_df = genes_df.drop("region", axis=1)
                genes_df = genes_df.drop("annot_genbank", axis=1)

                regions_df = self.crass_genes_df[['protein', 'region', 'annot_genbank']]
                genes_df = genes_df.merge(regions_df,
                                          left_on=genes_df.gene,
                                          right_on=regions_df.protein,
                                          how="inner")[['gene', 'gene_fam', 'annot_yutin_hmm', 'region', 'annot_genbank']]

            genes_df.loc[genes_df.annot_genbank.isnull(), "annot_genbank"] = "unknown"

            # add annotation from pvogs, we can add columns because MakeGenePlots
            merge_df = genes_df.merge(self.pvog_hmm_df,
                                      left_on=genes_df.gene,
                                      right_on=self.pvog_hmm_df.gene,
                                      how="left").drop(["key_0", "gene_y"], axis=1)
            merge_df.rename(columns={'gene_x': 'gene'}, inplace=True)

            merge_df = self.get_annotation_from_multiple_sources(merge_df)

            merge_df.to_csv(path_or_buf=gene_list_name, sep='\t', index=False)
            list_dfs.append(merge_df)

        concat_df = pd.concat(list_dfs)

        out_table_name2 = self.genome_dir + self.dir_sep + "out_gene_annotations_plus_pvog.txt"

        concat_df.to_csv(path_or_buf=out_table_name2, sep='\t', index=False)

    def get_annotation_from_multiple_sources(self, data):

        data.loc[data.annot_genbank == "hypothetical protein", "annot_genbank"] = "unknown"

        data["annotation"] = data.apply(self.derived_annotation, axis=1)

        return data

    @staticmethod
    def derived_annotation(row):
        annotation = row.annot_genbank
        source = "G"

        if annotation == "unknown":
            if not pd.isna(row.annot_yutin_hmm):
                annotation = row.annot_yutin_hmm
                source = "Y"

        if annotation == "Uncharacterized protein" or annotation == "unknown":
            if not pd.isna(row.annot_pvog):
                annotation = row.annot_pvog
                source = "P"

        return source + "-" + annotation


def annotate(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-ad", "--annotation_dir", dest="annotation_dir",
                        help="annotation directory with Yutin info", metavar="[annotation_dir]", required=True)

    parser.add_argument("-gd", "--genome_dir", dest="genome_dir",
                        help="genome directory with Guerin ref genomes", metavar="[genome_dir}", required=True)

    args = parser.parse_args(args_in)

    anno = AnnotateCrassGenomes(args.genome_dir, args.annotation_dir)

    anno.annotate()


# if __name__ == "__main__":
#     annotate(sys.argv[1:])

annotation_dir = r"D:\17 Dutihl Lab\source\phages_pycharm\input_files"
genome_dir = r"D:\17 Dutihl Lab\_tools\hmmsearch"

annotate(["-ad", annotation_dir, "-gd", genome_dir])

