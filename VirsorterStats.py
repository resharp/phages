import pandas as pd
import os
import re

class VirsorterStats:
    phage_signal_filename = ""
    all_affi_filename = ""

    phage_data = None
    all_affi_data = None
    phage_affi_data = None

    def __init__(self, mydir):

        self.phage_signal_filename = mydir + "VIRSorter_global-phage-signal.csv"
        self.all_affi_filename = mydir + "Metric_files\\VIRSorter_affi-contigs.tab"

        if not os.path.isfile(self.all_affi_filename):
            #switch to linux style
            self.all_affi_filename = mydir + "Metric_files/VIRSorter_affi-contigs.tab"

        self.read_phage_data(self.phage_signal_filename)
        self.read_affi_data(self.all_affi_filename)

    def read_phage_data(self, phage_signal_filename):

        f = open(phage_signal_filename, 'r')
        lines = f.readlines()
        header = lines[1].replace("## ", "")
        lines_no_comment = [n for n in lines if not n.startswith('##')]
        f.close()

        filename2 = phage_signal_filename + ".csv"
        f2 = open(filename2, 'w')
        f2.write(header)
        f2.writelines(lines_no_comment)
        f2.close()

        self.phage_data= pd.read_csv(filename2, delimiter=",")

        for i, row in self.phage_data.iterrows():
            fragment = self.phage_data.at[i, 'Fragment']

            if fragment.find("-gene_") > 0:
                phage_gene_start = re.sub("(.*)\-gene_([0-9]*)-gene_([0-9]*)", "\\2", fragment)
                phage_gene_end = re.sub("(.*)\-gene_([0-9]*)-gene_([0-9]*)", "\\3", fragment)
            else:
                phage_gene_start = 1
                phage_gene_end = int(self.phage_data.at[i, 'Nb genes'])

            self.phage_data.at[i, 'phage_gene_start'] = int(phage_gene_start)
            self.phage_data.at[i, 'phage_gene_end'] = int(phage_gene_end)

        # TODO: remove this copy if no longer necessary
        # self.phage_data.to_csv(filename2, sep=';', index=False)


    #TODO: change types of new fiels to int for better readability
    def read_affi_data(self, all_affi_filename):

        f = open(all_affi_filename, 'r')
        lines = f.readlines()
        lines_no_comment = [n for n in lines if not n.startswith('>')]
        f.close()

        filename2 = all_affi_filename + ".csv"
        f2 = open(filename2, 'w')
        f2.write("gene_id|nr1|nr2|nr3|nr4|gene_name|nr_6|nr_7|nr_8|short_annotation|nr_10|nr_11")
        f2.writelines(lines_no_comment)
        f2.close()

        self.all_affi_data = pd.read_csv(filename2, delimiter="|")
        self.all_affi_data = self.all_affi_data[self.all_affi_data.gene_name != ""]
        self.all_affi_data = self.all_affi_data[self.all_affi_data.gene_name != "-"]

        list_contigs = list(self.phage_data.Contig_id)

        #add derived field Contig_id for joining
        self.all_affi_data['Contig_id'] = self.all_affi_data.gene_id.str.replace("\-gene_[0-9]*", "")

        #now join all_affi_data and phage_data on Contig_id
        self.phage_affi_data = self.all_affi_data.merge(self.phage_data
                                             , left_on=self.all_affi_data.Contig_id
                                             , right_on=self.phage_data.Contig_id
                                             , how='inner')
        self.phage_affi_data = self.phage_affi_data[[
            'key_0'
            , 'gene_id'
            , 'gene_name'
            , 'short_annotation'
            , 'phage_gene_end'
            , 'phage_gene_start'
            , 'Nb genes'
        ]]
        self.phage_affi_data.rename(columns={'key_0':'Contig_id'}, inplace=True)

        for i, row in self.phage_affi_data.iterrows():
            gene_id = self.phage_affi_data.at[i, 'gene_id']
            gene_nr = self.gene_to_nr(gene_id)
            self.phage_affi_data.at[i, 'gene_nr'] = int(gene_nr)

        #now we only need the genes that are in the phage fragments
        self.phage_affi_data = self.phage_affi_data[
            (self.phage_affi_data.gene_nr >= self.phage_affi_data.phage_gene_start)
            & (self.phage_affi_data.gene_nr <= self.phage_affi_data.phage_gene_end)
        ]

        #TODO: remove this copy if no longer necessary
        # self.phage_affi_data.to_csv(filename2, sep=';', index=False)


    def category_counts(self):

        vcs = self.phage_data.Category.value_counts()

        vcs_dic = {}

        for i in vcs.index:
            vcs_dic[str(i)] = vcs[i]

        return vcs_dic

    def nr_of_phages(self):

        nr_of_phages = len(self.phage_data)

        return nr_of_phages

    def nr_of_complete_contigs(self):
        # -gene_298-gene_397

        #if the contig contains gene_ it seems to be a prophage

        #gene location in Fragment
        self.phage_data['Gene-loc'] = self.phage_data.Fragment.str.find("gene_")
        self.phage_data['Complete'] = self.phage_data['Gene-loc'] == -1

        # print("Non-complete ones")
        # print(self.phage_data[self.phage_data.Complete == False])

        return len(self.phage_data[self.phage_data.Complete])

    def all_gene_counts(self):
        return self.all_affi_data.gene_name.value_counts()

    def gene_to_contig(self, gene_id):
        contig_id = re.sub("\-gene_[0-9]*", "", gene_id)
        return contig_id

    def gene_to_nr(self, gene_id):
        contig_id = re.sub("(.*)\-gene_([0-9]*)", "\\2", gene_id)
        return contig_id

    def contig_gene_counts(self):
        contig_gene_counts = self.phage_affi_data.groupby(['Contig_id', 'gene_name']).size().reset_index()\
            .rename(columns={0: 'count'}).sort_values(['count'], ascending=[0])
        # df1.groupby(['A', 'B']).size().reset_index().rename(columns={0: 'count'})

        return contig_gene_counts


# self.phage_affi_data[self.phage_affi_data.gene_name.str.find("Phage_cluster") == 0].gene_name.value_counts()

