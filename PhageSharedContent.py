import pandas as pd
import os
import logging

# we would like to generate a list of [PH + PC] from [PH + IP]
# based on clustering of IPs (individual proteins) into PCs (protein clusters)

# then we would like to generate the shared PC content of all the combinations of phages
# in order to feed this into mcl
# the weight on each edge (connection between two phages) will be the number of shared genes

if os.name == "nt":
    mydir = r"D:\17 Dutihl Lab\_tools\mcl\1000"
else:
    mydir = "/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000"

extension = "mcl_75.I25"

class PhageSharedContent:

    phage_ip_table_name = "" #the file with the individual proteins per phage
    ip_pc_table_name = ""    #specific mcl result (clustering of IPs into PCs)
    phage_pc_table_name = ""
    shared_gene_content_name = ""
    pc_table_name = ""

    phage_ip_table = None
    ip_pc_table = None
    phage_pc_table = None
    pc_df = None

    def __init__(self, mydir, extension):

        logging.basicConfig(filename='PhageSharedContent.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.DEBUG)

        if os.name == 'nt':
            dir_sep = "\\"
        else:
            dir_sep = "/"

        #TODO: move file names to beginning of python script
        self.phage_ip_table_name = mydir + dir_sep + "phage_ip_table_short.txt"

        self.ip_pc_table_name = mydir + dir_sep + "ip_pc_table." + extension + ".filter_100"
        self.pc_table_name = mydir + dir_sep + "pc_table." + extension + ".filter_100"

        #output:
        self.phage_pc_table_name = mydir + dir_sep + "phage_pc_table." + extension + ".txt"
        self.shared_gene_content_name = mydir + dir_sep + "shared_gene_content." + extension + ".txt"


    def print_file_names(self):

        print("Start analyzing shared gene content between phages")
        print(self.phage_ip_table_name)
        print(self.ip_pc_table_name)
        print("Results will be put in:")
        print(self.shared_gene_content_name)

    def read_files(self):

        logging.debug("start reading tables")
        self.read_phage_ip_table()
        self.read_ip_pc_table()
        self.read_pc_table()
        logging.debug("finished reading tables")

    def read_phage_ip_table(self):
        self.phage_ip_table = pd.read_csv(self.phage_ip_table_name, delimiter=" ",
                                         header=None
                                         , usecols=[0, 1]
                                         , names=['phage_id', 'ip_id'])

    def read_ip_pc_table(self):
        self.ip_pc_table = pd.read_csv(self.ip_pc_table_name
                                         , delim_whitespace=True
                                         , header=None
                                         , usecols=[0, 1]
                                         , names=['ip_id', 'pc_id'])

    def read_pc_table(self):
        self.pc_df = pd.read_csv(self.pc_table_name, delimiter=" ",
                                 header=None
                                 , usecols=[0]
                                 , names=['pc_id'])

    def merge(self):
        self.phage_pc_table = self.phage_ip_table.merge(self.ip_pc_table
                                             , left_on=self.phage_ip_table.ip_id
                                             , right_on=self.ip_pc_table.ip_id
                                             , how='inner')

        self.phage_pc_table = self.phage_pc_table[['phage_id', 'pc_id']]

        self.phage_pc_table.to_csv(path_or_buf=self.phage_pc_table_name, index=False)

    def calc_shared_gene_content(self):

        # Dictionary of shared PCs between phages
        shared_pc_content = {}
        # shared_gene_content_dic[("PH1", "PH2")] = 0
        # shared_gene_content_dic[("PH2", "PH9")] = 3

        #now we can calculate the shared gene content
        self.phage_pc_table = self.phage_pc_table.sort_values('pc_id')

        # self.phage_pc_table.groupby(['pc_id']).size().reset_index().rename(columns={0: 'pc_count'}).sort_values(
        #     ['pc_count'], ascending=[0])

        #now loop through all clusters
        for index, row in self.pc_df.iterrows():
            pc_id = row["pc_id"]
            # with this pc_id make temporary dataframe filtered from phage content

            phage_df = self.phage_pc_table[ self.phage_pc_table.pc_id == pc_id]["phage_id"]

            # print(pc_id)
            # print(len(phage_df))

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
        f = open(self.shared_gene_content_name, "w+")
        for k, v in shared_pc_content.items():
            f.write("{} {} {}\n".format(k[0], k[1], v))
        f.close()

ph_content = PhageSharedContent(mydir, extension)

ph_content.print_file_names()

ph_content.read_files()

ph_content.merge()

ph_content.calc_shared_gene_content()


