import pandas as pd
import os
import logging
import math

# we would like to generate a list of [PH + PC] from [PH + IP]
# based on clustering of IPs (individual proteins) into PCs (protein clusters)

# then we would like to generate the shared PC content of all the combinations of phages
# in order to
# (A) feed this into mcl
#     the weight on each edge (connection between two phages) will be the number of shared genes
# (B) try to find a large group of phages with a minimum number of shared genes
#     therefore we would like to filter: what phages contain at least [min_nr_of_genes] with other phages
#

#TODO: Set skip_calculation to False in production situation
# * *  *   *    *      *        *
skip_shared_pc_calculation = True

#TODO: Set back to mcl_75.I25
# * *  *   *    *      *        *
#extension = "mcl_75.I25"
extension = "mcl_75.I20"

if os.name == "nt":
    mydir = r"D:\17 Dutihl Lab\_tools\mcl\1000"
else:
    mydir = "/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000"

min_nr_of_phages = 100  #how many genes shared between all of these phages?

class PhageSharedContent:

    min_nr_of_phages = 0

    phage_ip_table_name = "" #the file with the individual proteins per phage
    ip_pc_table_name = ""    #specific mcl result (clustering of IPs into PCs)
    phage_pc_table_name = ""
    phage_pc_table_name_filtered = ""

    shared_pc_content_name = ""
    shared_pc_measures_name = ""
    shared_phage_content_name = ""
    pc_table_name = ""
    phage_table_name = ""
    phage_table_name_filtered = ""

    phage_ip_table = None
    phage_ip_counts_df = None
    ip_pc_table = None
    phage_pc_table = None
    phage_pc_table_filtered = None

    pc_df = None
    shared_pcs_df = None
    shared_pcs_measures_df = None
    phage_df = None
    phage_df_filtered = None

    #for internal calculations
    max_nr_phages_possible = 0
    max_nr_of_pcs = 0

    def __init__(self, mydir, extension, min_nr_of_phages):

        logging.basicConfig(filename='PhageSharedContent.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.DEBUG)

        if os.name == 'nt':
            dir_sep = "\\"
        else:
            dir_sep = "/"

        #TODO: move file names to beginning of python script
        self.phage_ip_table_name = mydir + dir_sep + "phage_ip_table_short.txt"

        self.ip_pc_table_name = mydir + dir_sep + "ip_pc_table." + extension # + ".filter_100"

        # TODO: change
        self.pc_table_name = mydir + dir_sep + "pc_table." + extension # + ".filter_100"
        # self.pc_table_name = mydir + dir_sep + "pc_table." + extension + ".filter_100.for_testing"
        self.phage_table_name = mydir + dir_sep + "phage_table.txt"
        self.phage_table_name_filtered = mydir + dir_sep + "phage_table_sharing_{}.txt"

        #output:
        self.phage_pc_table_name = mydir + dir_sep + "phage_pc_table." + extension + ".txt"
        self.phage_pc_table_name_filtered = mydir + dir_sep + "phage_pc_table.sharing_{}." + extension + ".txt"

        self.shared_pc_content_name = mydir + dir_sep + "shared_gene_content." + extension + ".txt"
        self.shared_pc_measures_name = mydir + dir_sep + "shared_pc_measures." + extension + ".txt"

        self.shared_phage_content_name = mydir + dir_sep + "shared_phage_content." + extension + ".txt"


    def print_file_names(self):

        print("Start analyzing shared gene content between phages")
        print(self.phage_ip_table_name)
        print(self.ip_pc_table_name)
        print("Results will be put in:")
        print(self.shared_pc_content_name)
        print(self.shared_phage_content_name)

    def read_files(self):

        logging.debug("start reading tables")
        self.read_phage_ip_table()
        self.read_ip_pc_table()
        self.read_pc_table()
        self.read_phage_table()
        logging.debug("finished reading tables")

    def read_phage_ip_table(self):
        self.phage_ip_table = pd.read_csv(self.phage_ip_table_name, delimiter=" "
                                         , header=None
                                         , usecols=[0, 1]
                                         , names=['phage_id', 'ip_id'])

    def read_ip_pc_table(self):
        self.ip_pc_table = pd.read_csv(self.ip_pc_table_name
                                         , delim_whitespace=True
                                         , header=None
                                         , usecols=[0, 1]
                                         , names=['ip_id', 'pc_id'])

    def read_pc_table(self):
        self.pc_df = pd.read_csv(self.pc_table_name, delimiter=" "
                                 , header=None
                                 , usecols=[0]
                                 , names=['pc_id'])

    def read_phage_table(self):
        self.phage_df = pd.read_csv(self.phage_table_name
                                 , delim_whitespace=True
                                 , header=None
                                 , usecols=[0,1]
                                 , names=['ph_id', 'ph_name'])

    def merge(self):
        self.phage_pc_table = self.phage_ip_table.merge(self.ip_pc_table
                                             , left_on=self.phage_ip_table.ip_id
                                             , right_on=self.ip_pc_table.ip_id
                                             , how='inner')

        self.phage_pc_table = self.phage_pc_table[['phage_id', 'pc_id']]

        self.phage_pc_table.to_csv(path_or_buf=self.phage_pc_table_name, index=False)

    def calc_shared_pc_content(self):

        # Dictionary of shared PCs between phages
        shared_pc_content = {}
        # shared_gene_content_dic[("PH1", "PH2")] = 0
        # shared_gene_content_dic[("PH2", "PH9")] = 3

        #now we can calculate the shared gene content
        self.phage_pc_table = self.phage_pc_table.sort_values('pc_id')

        #now loop through all clusters
        for index, row in self.pc_df.iterrows():
            pc_id = row["pc_id"]
            # with this pc_id make temporary dataframe filtered from phage content

            #TODO either
            #  Fix bug to ignore double entries
            # OR ignore these PCs (or phages?)
            phage_df = self.phage_pc_table[ self.phage_pc_table.pc_id == pc_id]["phage_id"].unique()

            index = pd.MultiIndex.from_product([phage_df, phage_df], names = ["phage_1", "phage_2"])

            phage_phage_df = pd.DataFrame(index = index).reset_index()
            print('calculating for protein cluster: ' + str(pc_id) + ' shared by ' + str(len(phage_df)) + ' phages.')
            print(len(phage_phage_df))

            for index, row in phage_phage_df.iterrows():
                phage_1 = row["phage_1"]
                phage_2 = row["phage_2"]

                if phage_1 != phage_2:
                    key = (phage_1, phage_2)
                    if key in shared_pc_content:
                        shared_pc_content[key] += 1
                    else:
                        shared_pc_content[key] = 1

        #now save the shared gene content to file
        rows_list = []

        for k, v in shared_pc_content.items():
            rows_list.append([k[0], k[1], v])
        self.shared_pcs_df = pd.DataFrame(rows_list, columns=["phage_1", "phage_2", "count"])

        self.save_shared_pc_content()

    def save_shared_pc_content(self):
        self.shared_pcs_df.to_csv(path_or_buf=self.shared_pc_content_name
                                  , index=False, header=None)

    def print_statistics_of_gene_content(self):

        nr_phages = len(self.phage_df)
        log_line = "nr of phages: {}. Request for group of at least {} phages".\
            format(nr_phages, self.min_nr_of_phages)
        logging.info(log_line)
        print(log_line)

        self.max_nr_phages_possible = 0

        # we would like to filter: what phages contain at least [min_nr_of_genes] with other phages
        # and consist of a group of min_nr_of_phages

        print_cutoff = False

        #TODO: change 1000 to length of largest phage (in # of pcs)
        for min_nr_of_pcs in range(1,1000):

            #note: divide by two because of double (reversed) entries
            nr_phages_combinations = int(len(self.shared_pcs_df[self.shared_pcs_df["count"] > min_nr_of_pcs - 1] ) / 2)
            max_nr_phages_possible = int(math.sqrt(nr_phages_combinations))


            nr_phages = len(self.shared_pcs_df[self.shared_pcs_df["count"] > min_nr_of_pcs - 1].groupby(["phage_1"]).count())

            if max_nr_phages_possible < self.min_nr_of_phages:
                if not print_cutoff:
                    logging.debug("---------------------")
                    logging.debug("maximum reached: {} phages with max nr of pcs {}.".format(
                        self.max_nr_phages_possible
                        , self.max_nr_of_pcs))
                    print_cutoff = True
            else:
                self.max_nr_phages_possible = max_nr_phages_possible
                self.max_nr_of_pcs = min_nr_of_pcs

            log_line = "# phage combinations sharing {} pcs or more: {}, consisting of {} phages, max phages sharing all genes: {} ".format(
                min_nr_of_pcs
            ,   nr_phages_combinations
            ,   nr_phages
            ,   max_nr_phages_possible)
            logging.debug(log_line)
            print(log_line)

            #TODO: what are the number of phages in this set?
            #TODO: the number of phage combinations should be larger than min_nr_of_phages^2

            if nr_phages_combinations == 0:
                break

    def read_shared_pc_content(self):
        self.shared_pcs_df = pd.read_csv(self.shared_pc_content_name
                                         , delimiter=","
                                         , header=None
                                         , usecols=[0,1,2]
                                         , names=["phage_1", "phage_2", "count"])

    def determine_largest_cluster_sharing_n_pcs(self, nr_of_pc_shared):

        # we want to have a list of phages containing one gene
        # for this we can just take all
        # phages from
        # self.shared_pcs_df if it is df
        if nr_of_pc_shared == 1:
            shared_pcs_df_filtered = self.shared_pcs_df
        else:
            shared_pcs_df_filtered = self.shared_pcs_df[self.shared_pcs_df["count"] > nr_of_pc_shared - 1 ]

        self.phage_df_filtered = shared_pcs_df_filtered.groupby(["phage_1"]).count().add_suffix('_Count').reset_index()[['phage_1']]

        self.phage_pc_table_filtered = self.phage_pc_table.merge(self.phage_df_filtered
                                             , left_on=self.phage_pc_table.phage_id
                                             , right_on=self.phage_df_filtered["phage_1"]
                                             , how='inner')[['phage_id','pc_id']]
        # join with phage table (for also having the names)
        print("number of phages")
        print(len(self.phage_df_filtered))

        print("length of phage pc table")
        print(len(self.phage_pc_table_filtered))

        print("length of phage pc table after filtering")
        print(len(self.phage_pc_table_filtered))

        self.phage_df_filtered = self.phage_df_filtered.merge(self.phage_df
                                            , left_on = self.phage_df_filtered.phage_1
                                            , right_on =self.phage_df.ph_id
                                            , how='inner') [['ph_id', 'ph_name']]

        self.phage_df_filtered.to_csv(path_or_buf=self.phage_table_name_filtered.format(nr_of_pc_shared)
                                  , index=False, header=None)
        self.phage_pc_table_filtered.to_csv(path_or_buf=self.phage_pc_table_name_filtered.format(nr_of_pc_shared)
                                  , index=False, header=None)

    def filter_on_thresholds(self):

        max_nr = self.max_nr_of_pcs

        # determine list of phages based on cutoff
        shared_pcs_df_filtered = self.shared_pcs_df[self.shared_pcs_df["count"] > max_nr - 1 ]

        #get the phages out
        self.phage_df_filtered_= shared_pcs_df_filtered.groupby(["phage_1"]).count().add_suffix('_Count').reset_index()[['phage_1']]

        print("number of phages")
        print(len(self.phage_df_filtered))

        print("length of phage pc table")
        print(len(self.phage_pc_table))

        self.phage_pc_table_filtered = self.phage_pc_table.merge(self.phage_df_filtered
                                             , left_on=self.phage_pc_table.phage_id
                                             , right_on=self.phage_df_filtered["phage_1"]
                                             , how='inner')[['phage_id','pc_id']]

        print("length of phage pc table after filter")
        print(len(self.phage_pc_table_filtered))

    def calc_shared_pc_measures(self):

        #now calculate phage_ip_counts from phage_ip_table

        self.phage_ip_counts_df = self.phage_ip_table.groupby(["phage_id"]).count().add_suffix('_Count').reset_index()
        #count column will be ip_id_Count

        #optional: write phage_ip_counts

        #join shared_pcs_df on phage_1 column
        self.shared_pcs_measures_df = self.shared_pcs_df.merge(self.phage_ip_counts_df
                                 , left_on=self.shared_pcs_df.phage_1
                                 , right_on=self.phage_ip_counts_df.phage_id
                                 , how='inner')
        self.shared_pcs_measures_df.rename(columns={'ip_id_Count':'phage_1_ip_count'}, inplace=True)
        self.shared_pcs_measures_df = self.shared_pcs_measures_df[['phage_1','phage_2','count','phage_1_ip_count']]

        #join shared_pcs_df on phage_2 column

        self.shared_pcs_measures_df = self.shared_pcs_measures_df.merge(self.phage_ip_counts_df
                                 , left_on=self.shared_pcs_measures_df.phage_2
                                 , right_on=self.phage_ip_counts_df.phage_id
                                 , how='inner')
        self.shared_pcs_measures_df.rename(columns={'ip_id_Count':'phage_2_ip_count'}, inplace=True)

        self.shared_pcs_measures_df = self.shared_pcs_measures_df[['phage_1','phage_2','count','phage_1_ip_count','phage_2_ip_count']]

        #add extra column for computation of Jaccard similarity score (use apply)
        self.shared_pcs_measures_df['jaccard_index'] = self.shared_pcs_measures_df['count'] / \
                                                        (self.shared_pcs_measures_df['phage_1_ip_count'] +\
                                                        self.shared_pcs_measures_df['phage_2_ip_count'] -\
                                                        self.shared_pcs_measures_df['count'])

        self.shared_pcs_measures_df = self.shared_pcs_measures_df.sort_values(['jaccard_index'], ascending=[0])

        #write shared_pcs_measures
        self.shared_pcs_measures_df.to_csv(path_or_buf=self.shared_pc_measures_name, index=False)

ph_content = PhageSharedContent(mydir, extension, min_nr_of_phages)

ph_content.print_file_names()

ph_content.read_files()

ph_content.merge()

if not skip_shared_pc_calculation:
    # this one fills self.shared_pcs_df, if already done before you can skip this while testing
    ph_content.calc_shared_pc_content()
else:
    ph_content.read_shared_pc_content()

ph_content.print_statistics_of_gene_content()

for i in range(1,2):
    ph_content.determine_largest_cluster_sharing_n_pcs(i)

ph_content.calc_shared_pc_measures()

#ph_content.filter_on_thresholds()

#Note: it only makes sense to do an iteration if you filtered both on
# min_nr_of_pcs AND on
# .. min_nr_of_phages

# ph_content.calc_shared_gene_content()
# ph_content.print_statistics_of_gene_content()
