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
    mydir = "/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000"

extension = "mcl_75.I25" # extension for mcl ip clustering
extension2 = "I25" # extension for mcl phage clustering

class PhageContentPlots:

    extension = ""
    mydir = ""

    ip_pc_table_name = ""    #specific mcl result (clustering of IPs into PCs)
    # ip_pc_table.mcl_75.I20

    phc_ph_table_name = ""   #clustering of phages in phage clusters
    # TODO: reformat in other format (we need a coupling between PHC en PH, it is not in that format yet)
    # out.shared_gene_content.mcl_75.I20.txt.I25

    ip_pc_table = None

    def __init__(self, mydir, extension):

        self.mydir = mydir
        self.extension = extension

        if os.name == 'nt':
            dir_sep = "\\"
        else:
            dir_sep = "/"

        self.ip_pc_table_name = mydir + dir_sep + "ip_pc_table." + extension
        # self.phc_ph_table_name = mydir + dir_sep + "[TODO]." + extension

    def read_ip_pc_table(self):
        self.ip_pc_table = pd.read_csv(self.ip_pc_table_name
                                         , delim_whitespace=True
                                         , header=None
                                         , usecols=[0, 1]
                                         , names=['ip_id', 'pc_id'])

    def make_pc_distribution(self):
        self.read_ip_pc_table()

        #here we should prepare the data
        data = self.ip_pc_table.groupby(["pc_id"]).count().add_suffix('_Count').reset_index()
        max_size = data.iloc[0].values[1]
        data = data[["ip_id_Count"]]

        data = data["ip_id_Count"].values.tolist()

        min_ips_in_pcs = 100
        max_ips_in_pcs = max_size

        self.make_pc_dis_with_cutoff(data, 2, max_ips_in_pcs)
        self.make_pc_dis_with_cutoff(data, min_ips_in_pcs, max_ips_in_pcs)

    def make_pc_dis_with_cutoff(self, data, min_ips_in_pcs, max_ips_in_pcs):

        plt.clf()
        plt.hist(data, bins=range(min_ips_in_pcs, max_ips_in_pcs))
        plt.xlabel("nr. of proteins in cluster")
        plt.ylabel("frequency")
        plt.title("Occurence of protein cluster sizes for {}\n "
                  "minimum ips in pc set at: {}".format(self.extension, min_ips_in_pcs))
        # plt.show()
        if os.name == 'nt':
            dir_sep = "\\"
        else:
            dir_sep = "/"
        figure_name = "{}{}ph_plots.pc_distribution.min_{}.{}.png".format(self.mydir, dir_sep, min_ips_in_pcs, self.extension)
        plt.savefig(figure_name)

ph_plots = PhageContentPlots(mydir, extension)

ph_plots.make_pc_distribution()
