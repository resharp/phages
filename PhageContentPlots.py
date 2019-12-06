import pandas as pd
import os
import matplotlib.pyplot as plt

# we would like to show distributions of
# 1. # ips (proteins) in pcs (protein clusters)
# 2. # ips (proteins) in pccs (protein cluster clusters)
# 3. # viruses in virus clusters
#
#
# only read in the files that you need
# I might also add a test class for this
if os.name == "nt":
    mydir = r"D:\17 Dutihl Lab\_tools\mcl\1000"
else:
    mydir = "/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_all"

extension = "mcl_75.I25" # extension for mcl ip clustering
extension2 = "I25" # extension for mcl phage clustering

class PhageContentPlots:

    extension = ""
    mydir = ""
    dir_sep = ""

    ip_pc_table_name = ""    #specific mcl result (clustering of IPs into PCs)
    # ip_pc_table.mcl_75.I20

    phc_ph_table_name = ""   #clustering of phages in phage clusters
    # TODO: reformat in other format (we need a coupling between PHC en PH, it is not in that format yet)
    # out.shared_gene_content.mcl_75.I20.txt.I25

    phage_ip_table_name = "" #the file with the individual proteins per phage

    shared_pc_measures_name = ""

    ip_pc_table = None
    phage_ip_table = None
    shared_pcs_measures_df = None

    def __init__(self, mydir, extension):

        self.mydir = mydir
        self.extension = extension

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        self.ip_pc_table_name = mydir + self.dir_sep + "ip_pc_table." + extension
        self.phage_ip_table_name = mydir + self.dir_sep + "phage_ip_table_short.txt"
        # self.phc_ph_table_name = mydir + self.dir_sep + "[TODO]." + extension
        self.shared_pc_measures_name = mydir + self.dir_sep + "shared_pc_measures." + extension + ".txt"


    def read_ip_pc_table(self):
        self.ip_pc_table = pd.read_csv(self.ip_pc_table_name
                                         , delim_whitespace=True
                                         , header=None
                                         , usecols=[0, 1]
                                         , names=['ip_id', 'pc_id'])

    def read_phage_ip_table(self):
        self.phage_ip_table = pd.read_csv(self.phage_ip_table_name
                                          ,delim_whitespace=True
                                          , header=None
                                          , usecols=[0,1]
                                          , names=['phage_id', 'ip_id'])
    def read_shared_pc_measures(self):
        self.shared_pcs_measures_df= pd.read_csv(self.shared_pc_measures_name
                                        ,delimiter=","
                                        ,usecols=[0,1,2,3,4,5]
                                        ,dtype={"jaccard_index" : "float64"})

        dummy = "true"

    def make_pc_distribution(self):
        self.read_ip_pc_table()

        data = self.ip_pc_table.groupby(["pc_id"]).count().add_suffix('_Count').reset_index()
        max_size = data.iloc[0].values[1]
        data = data[["ip_id_Count"]]

        data = data["ip_id_Count"].values.tolist()

        max_ips_in_pcs = max_size

        self.make_pc_dis_with_cutoff(data, 2, max_ips_in_pcs + 100)

    def make_pc_dis_with_cutoff(self, data, min_ips_in_pcs, max_ips_in_pcs):

        plt.clf()
        plt.hist(data, bins=range(min_ips_in_pcs, max_ips_in_pcs, 1), log=True)
        plt.xlabel("nr. of proteins in cluster")
        plt.ylabel("log10 of frequency")
        plt.title("Occurence of protein cluster sizes for {}\n "
                  "minimum ips in pc set at: {}".format(self.extension, min_ips_in_pcs))

        figure_name = "{}{}ph_plots.pc_distribution.min_{}.{}.pdf".format(self.mydir, self.dir_sep, min_ips_in_pcs, self.extension)
        plt.savefig(figure_name)

    def make_phage_distribution(self):

        self.read_phage_ip_table()

        data = self.phage_ip_table.groupby(["phage_id"]).count().add_suffix('_Count').reset_index()
        data = data[["ip_id_Count"]]

        data = data["ip_id_Count"].values.tolist()

        max_size = max(data)
        plt.clf()

        plt.hist(data, bins=range(1, max_size, 1), log=True)

        plt.xlabel("nr. of ips (proteins) in (pro)phage")
        plt.ylabel("log10 of frequency")
        plt.title("Distribution of number of proteins in predicted phages.")

        figure_name = "{}{}ph_plots.phage_distribution.pdf".format(self.mydir, self.dir_sep)
        plt.savefig(figure_name)
        # plt.show()

    def make_jaccard_distribution(self):

        ph_plots.read_shared_pc_measures()
        data = self.shared_pcs_measures_df

        data = data['jaccard_index'].values.tolist()

        plt.clf()
        plt.xlabel("Jaccard similarity coefficient")
        plt.ylabel("log10 scale of frequency")
        plt.title("Distribution of pairwise similarities of phages for {}".format(self.extension))

        plt.hist(data, log=True, bins=[ i/200 for i in range(1,201)])

        figure_name = "{}{}ph_plots.jaccard_distribution.{}.pdf".format(self.mydir, self.dir_sep, self.extension)
        plt.savefig(figure_name)


ph_plots = PhageContentPlots(mydir, extension)

ph_plots.make_pc_distribution()

#TODO: also make a distribution of the number of pcs that are in the phages (because not all ips are in cps)
ph_plots.make_phage_distribution()

ph_plots.make_jaccard_distribution()
