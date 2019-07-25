
import pandas as pd
import os
import re
import logging
import io

class VirsorterStats:
    phage_signal_filename = ""
    all_affi_filename = ""

    phage_data = None
    all_affi_data = None
    phage_affi_data = None
    gene_counts = None

    def __init__(self, mydir):

        self.phage_signal_filename = mydir + "VIRSorter_global-phage-signal.csv"
        self.all_affi_filename = mydir + "Metric_files\\VIRSorter_affi-contigs.tab"

        if not os.path.isfile(self.all_affi_filename):
            #switch to linux style
            self.all_affi_filename = mydir + "Metric_files/VIRSorter_affi-contigs.tab"

        logging.basicConfig(filename='VirsorterStats.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.DEBUG)
        logging.info('Start VirsorterStats')

        self.read_phage_data(self.phage_signal_filename)

        self.read_affi_data(self.all_affi_filename)

        self.calc_gene_counts(mydir)

    def read_phage_data(self, phage_signal_filename):

        logging.debug("start opening file with phage data")
        f = open(phage_signal_filename, 'r')

        logging.debug("start reading all lines")
        lines = f.readlines()
        header = lines[1].replace("## ", "")
        logging.debug("start filtering out comments in all lines with list comprehension")
        lines_no_comment = [n for n in lines if not n.startswith('##')]
        f.close()

        filename2 = phage_signal_filename + ".csv"
        f2 = open(filename2, 'w')
        f2.write(header)
        logging.debug("start writing all lines to intermediary file")
        f2.writelines(lines_no_comment)
        logging.debug("end writing all lines to intermediary file")
        f2.close()

        logging.debug("start read phage data from csv")
        self.phage_data= pd.read_csv(filename2, delimiter=",")
        logging.debug("end read phage data from csv")

        os.remove(filename2)

        logging.debug("start looping and doing string operations on every single row")
        for i, row in self.phage_data.iterrows():
            fragment = self.phage_data.at[i, 'Fragment']

            #TODO Improve performance Can we do this in one go instead of three?
            if fragment.find("-gene_") > 0:
                phage_gene_start = re.sub("(.*)\-gene_([0-9]*)-gene_([0-9]*)", "\\2", fragment)
                phage_gene_end = re.sub("(.*)\-gene_([0-9]*)-gene_([0-9]*)", "\\3", fragment)
            else:
                phage_gene_start = 1
                phage_gene_end = int(self.phage_data.at[i, 'Nb genes'])

            self.phage_data.at[i, 'phage_gene_start'] = int(phage_gene_start)
            self.phage_data.at[i, 'phage_gene_end'] = int(phage_gene_end)
        logging.debug("end looping and doing string operations on every single row")

        # TODO: remove this copy if no longer necessary
        # self.phage_data.to_csv(filename2, sep=';', index=False)


    #TODO: change types of new fiels to int for better readability
    def read_affi_data(self, all_affi_filename):

        logging.debug("start open affi data")
        f = open(all_affi_filename, 'r')
        logging.debug("start read all lines of affi data in one go")
        lines = f.readlines()

        ##TODO: Improve performance, can we read the lines in chunks and not load everything in memory
        logging.debug("start filtering all lines on comments")
        #TODO: do this with grep it is way faster
        lines_no_comment = [n for n in lines if not n[0] == '>'   #quicker than startswith
                            and not n.__contains__('|-|-|-|-|-|') #this line quickly leaves out empty annotations
                            ]
        #lines_no_comment = [n for n in lines if not n[0] == '>']
        f.close()

        filename2 = all_affi_filename + ".csv"
        f2 = open(filename2, 'w')

        logging.debug("start writing all lines of affi data to intermediary file")
        f2.writelines(lines_no_comment)
        filename2 = all_affi_filename + ".csv"
        f2.close()

        logging.debug("start read affi data from csv")

        #TODO: read in chunks to improve performance and memory usage
        self.all_affi_data = pd.read_csv(filename2, delimiter="|",
                                         header=None
                                         , usecols=[0,5,9]
                                         , names=['gene_id','gene_name','short_annotation'])

        os.remove(filename2)

        logging.debug("start filtering out empty gene names")
        self.all_affi_data = self.all_affi_data[(self.all_affi_data.gene_name != "") & (self.all_affi_data.gene_name != "-")]
        logging.debug("end filtering out empty gene names")

        #add derived field Contig_id for joining

        logging.debug("start calculating Contig_id")
        self.all_affi_data['Contig_id'] = self.all_affi_data.gene_id.str.replace("\-gene_[0-9]*", "")

        #now join all_affi_data and phage_data on Contig_id
        logging.debug("start merging with phage data")
        #TODO: would it be possible to limit the number of columns beforehand?
        self.phage_affi_data = self.all_affi_data.merge(self.phage_data
                                             , left_on=self.all_affi_data.Contig_id
                                             , right_on=self.phage_data.Contig_id
                                             , how='inner')

        del self.all_affi_data #remove largest in memory object asap

        logging.debug("start selecting just af few columns of the affi data")
        self.phage_affi_data = self.phage_affi_data[[
            'key_0'
            , 'gene_id'
            , 'gene_name'
            , 'short_annotation'
            , 'phage_gene_end'
            , 'phage_gene_start'
            , 'Nb genes'
        ]]
        logging.debug("start renaming column")
        self.phage_affi_data.rename(columns={'key_0':'Contig_id'}, inplace=True)

        logging.debug("start computing gene_nr with regex")
        self.phage_affi_data['gene_nr'] = self.phage_affi_data.gene_id.apply(self.gene_to_nr)

        #now we only need the genes that are in the phage fragments

        logging.debug("start filtering out genes that are not in the phage fragments")
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

    # def all_gene_counts(self):
    #     return self.all_affi_data.gene_name.value_counts()

    def gene_to_contig(self, gene_id):
        contig_id = re.sub("\-gene_[0-9]*", "", gene_id)
        return contig_id

    def gene_to_nr(self, gene_id):

        contig_id = re.sub("(.*)\-gene_([0-9]*)", "\\2", gene_id)

        return int(contig_id)

    def get_gene_counts_df(self):
        contig_gene_counts = self.phage_affi_data.groupby(['gene_name']).size().reset_index()\
            .rename(columns={0: 'gene_count'}).sort_values(['gene_count'], ascending=[0])
        # df1.groupby(['A', 'B']).size().reset_index().rename(columns={0: 'count'})

        return contig_gene_counts

    def calc_gene_counts(self, mydir):

        self.gene_counts = self.get_gene_counts_df()
        filename = mydir + "shared_gene_counts.csv"

        self.gene_counts = self.gene_counts[self.gene_counts.gene_count > 1]

        self.gene_counts = self.gene_counts.merge(self.phage_affi_data
                                             , left_on=self.gene_counts.gene_name
                                             , right_on=self.phage_affi_data.gene_name
                                             , how='inner')

        df = self.gene_counts



        self.gene_counts = pd.DataFrame({'count' : df.groupby( ['key_0', 'gene_count', 'short_annotation'] ).size()})\
            .reset_index() \
            .sort_values(['gene_count'], ascending=[0])[['key_0', 'gene_count', 'short_annotation']]

        self.gene_counts.rename(columns={'key_0':'gene_name'}, inplace=True)

        self.gene_counts.to_csv(path_or_buf=filename, index=False)





# self.phage_affi_data[self.phage_affi_data.gene_name.str.find("Phage_cluster") == 0].gene_name.value_counts()

