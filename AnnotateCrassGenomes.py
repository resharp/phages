import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import argparse

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

    dir_sep = ""

    ref_genes_df = None
    hmm_hits_df = None
    yutin_genes_df = None
    crass_genes_df = None

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

        pass

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
                                       , names=["gene", "gene_fam", "hmm_evalue", "hmm_score"
                                                , "start","end","AA_length"]
                                       , dtype={"start": "int32", "end": "int32", "AA_length": "int32"}
                                       )

        # print(self.hmm_hits_df.dtypes)
        print("number of hmm hits: {nr_hits}".format(nr_hits=len(self.hmm_hits_df)))

        self.yutin_genes_df = pd.read_csv(self.yutin_genes_name
                                          , sep="\t"
                                          , index_col=None
                                          , usecols=range(0,23)
                                          )
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

        merge_df = merge_df.merge(self.yutin_genes_df[["gene_fam","yutin_gene_nr","gene_annot"]],
                                  left_on=merge_df.gene_fam,
                                  right_on=self.yutin_genes_df.gene_fam,
                                  how="left").drop(["key_0"], axis=1)
        merge_df = merge_df.merge(self.crass_genes_df,
                                  left_on=merge_df.yutin_gene_nr,
                                  right_on=self.crass_genes_df.yutin_gene_nr,
                                  how="left").drop(["key_0"], axis=1)
        debug = "True"

        out_table_name = self.genome_dir+ self.dir_sep + "out_gene_annotations.txt"

        merge_df.to_csv(path_or_buf=out_table_name, sep='\t', index=False)


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

# TODO for testing, do not use in production
annotation_dir = r"D:\17 Dutihl Lab\crassphage"
genome_dir = r"D:\17 Dutihl Lab\_tools\hmmsearch"

annotate(["-ad", annotation_dir, "-gd", genome_dir])

