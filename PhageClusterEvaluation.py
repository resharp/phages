import pandas as pd
import os
import matplotlib.pyplot as plt
import igraph as ig

#TODO create directory, e.g. /mcl_75.I25.plots

# we would like to visually and automatically inspect
# the content of phage cluster
# - how large is the cluster?
# - which hosts do they come from?
# - what is the Jaccard distribution of the clusters? (do we see a structure? multiple peaks?)
# - what is the plain phage with the highest average (Jaccard) similarity to all other phages in the cluster?
#       this is not strictly necessary we might also do readmapping against all phages in the cluster
# - distribution of Jaccard index from the plain phage?
# - can we make an iGraph plot of each cluster to assess structure?
#   use e.g. Fruchterman-Reingold layout, it should use the weight of edges as strength
# - megablast the plain phage

# INPUT
# genome_phage_table_short.txt                                          (coupled PATRIC host id and phage predictions)
# genome_taxonomy.txt                                                   (host taxonomy from PATRIC)
# shared_pc_measures.mcl_75.I25.txt                                     (Jaccard similarities and nr. of ips_
# split.shared_pc_measures.mcl_75.I25.nr_pc_min_5.sharing_0_1.abc.I25   (phage clusters)

if os.name == "nt":
    mydir = r"D:\17 Dutihl Lab\_tools\mcl\all"
else:
    mydir = "/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000"

extension = "mcl_75.I25" # extension for mcl ip clustering (create subdir if it does not exist)
extension2 = "I25" # extension for mcl phage clustering

class PhageClusterEvaluation:

    extension = ""
    mydir = ""
    dir_sep = ""

    ip_pc_table_name = ""    #specific mcl result (clustering of IPs into PCs)
    # ip_pc_table.mcl_75.I20

    phc_ph_table_name = ""   #clustering of phages in phage clusters

    ph_distances_name = "" #the file with the pairwise Jaccard similarities (and phage sizes in nr of IPs)

    genome_ph_name = ""
    genome_taxonomy_name = ""

    ph_distances_df = None
    phc_ph_df = None
    genome_ph_df = None
    genome_taxonomy_df = None

    all_df = None

    def __init__(self, mydir, extension, extension2):

        self.mydir = mydir
        self.extension = extension
        self.extension2 = extension2

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        #self.ip_pc_table_name = mydir + self.dir_sep + "ip_pc_table." + extension
        self.phc_ph_table_name = mydir + self.dir_sep + "split.shared_pc_measures.mcl_75.I25.nr_pc_min_5.sharing_05.abc." + extension2
        self.ph_distances_name = mydir + self.dir_sep + "shared_pc_measures." + extension + ".nr_pc_min_5.sharing_05.abc"

        self.genome_ph_name = mydir + self.dir_sep + "genome_phage_table_short.txt"
        self.genome_taxonomy_name = mydir + self.dir_sep + "genome_taxonomy.txt"

    def read_files(self):

        self.ph_distances_df= pd.read_csv(self.ph_distances_name, delimiter=" "
                                        , skiprows=[0]
                                        , usecols=[0, 1, 2]
                                        , names=[ 'phage_1'
                                                , 'phage_2'
                                                # , 'count'
                                                # , 'phage_1_ip_count'
                                                # , 'phage_2_ip_count'
                                                , 'jaccard_index'])

        self.phc_ph_df= pd.read_csv(self.phc_ph_table_name, delimiter=" "
                                        , header=None
                                        , usecols=[0, 1]
                                        , names=[ 'phc_id'
                                                , 'phage_id'])

        self.genome_ph_df= pd.read_csv(self.genome_ph_name, delimiter=" "
                                        , header=None
                                        , usecols=[0, 1]
                                        , names=[ 'genome_id'
                                                , 'phage_id'])

        self.genome_taxonomy_df= pd.read_csv(self.genome_taxonomy_name
                                        , sep='\t'
                                        , header=None
                                        , usecols=[0, 1, 2, 3]
                                        , names=[ 'genome_id'
                                                , 'family'
                                                , 'genus'
                                                , 'species']
                                        , dtype={'genome_id' : 'object'})
        # self.genome_taxonomy_df.info()

    def determine_origin(self):
        # for each cluster:
        # - which hosts do they come from?

        self.phc_ph_df.merge

        merge_df = self.phc_ph_df.merge(self.genome_ph_df
                                     , left_on=self.phc_ph_df.phage_id
                                     , right_on=self.genome_ph_df.phage_id
                                     , how='inner')

        merge_df = merge_df[['phc_id', 'phage_id_x', 'genome_id']]

        merge_df.rename(columns={'phage_id_x': 'phage_id'}, inplace=True)

        merge_df = merge_df.merge(self.genome_taxonomy_df
                                  , left_on=merge_df.genome_id
                                  , right_on=self.genome_taxonomy_df.genome_id
                                  , how='inner')

        merge_df = merge_df[['phc_id','phage_id','genome_id_x','family','genus','species']]

        merge_df.rename(columns={'genome_id_x': 'genome_id'}, inplace=True)

        self.all_df = merge_df


    def make_jaccard_distribution(self):

        phcs = []

        for i in range(1,2):
        # for i in range(1,11):
            phcs.append("PHC_" + str(i))

        for phc in phcs:
            self.jaccard_dis_for_phage_cluster(phc)

    def jaccard_dis_for_phage_cluster(self, phc):
        data = self.phc_ph_df[self.phc_ph_df.phc_id == phc][['phage_id']]

        #join distances with the phages of one cluster
        data = data.merge(self.ph_distances_df
                          , left_on=data.phage_id
                          , right_on=self.ph_distances_df.phage_1
                          , how='inner')

        plt.clf()
        plot = data['jaccard_index'].plot(kind='density', bw_method=0.05)

        plt.xlabel("Jaccard similarity coefficient")
        plt.ylabel("density")
        plt.title("Pairwise similarities of phages for cluster {}, inflation {}".format(phc, self.extension2))

        figure_name = "{}{}{}{}ph_plots.jaccard_distribution_for_cluster.{}.density.{}.pdf" \
            .format(self.mydir, self.dir_sep, "mcl_75.I25.plots", self.dir_sep, phc, self.extension2)

        plot.get_figure().savefig(figure_name, format='pdf')

        data = data['jaccard_index'].values.tolist()
        plt.clf()
        plt.xlabel("Jaccard similarity coefficient")
        plt.ylabel("log10 scale of frequency")
        plt.title("Pairwise similarities of phages for cluster {}, inflation {}".format(phc, self.extension2))

        plt.hist(data, log=True, bins=[i / 100 for i in range(50, 101)])
        figure_name = "{}{}{}{}ph_plots.jaccard_distribution_for_cluster.{}.{}.pdf" \
            .format(self.mydir, self.dir_sep, "mcl_75.I25.plots", self.dir_sep, phc, self.extension2)
        # plt.show()
        plt.savefig(figure_name)

    def make_cluster_graphs(self):

        phcs = []
        # for i in range(1,2):
        for i in range(1,11):
            phcs.append("PHC_" + str(i))

        for phc in phcs:
            self.make_cluster_graph(phc)

    def make_cluster_graph(self, phc):

        #now collect all pairwise interactions for phc

        phages_df = self.phc_ph_df[self.phc_ph_df.phc_id == phc][['phage_id']]

        cluster_dist_df = phages_df.merge(self.ph_distances_df
                            ,left_on=phages_df.phage_id
                            ,right_on=self.ph_distances_df.phage_1
                            ,how='inner')[['phage_1', 'phage_2', 'jaccard_index']]

        #also filter on links towards phages in other clusters
        cluster_dist_df = phages_df.merge(cluster_dist_df
                            ,left_on=phages_df.phage_id
                            ,right_on=cluster_dist_df.phage_2
                            ,how='inner')[['phage_1', 'phage_2', 'jaccard_index']]

        g = ig.Graph()

        #first determine all vertices

        phages_list = phages_df['phage_id'].values.tolist()

        #TODO: you can also do this in one go, without looping
        for phage in phages_list:
            g.add_vertices(phage)

        for index, row in cluster_dist_df.iterrows():
            phage_1 = row.phage_1
            phage_2 = row.phage_2

            #check if the edge in the other direction does not exist!
            if g.get_eid(phage_2, phage_1, directed=False, error=False) == -1:
                g.add_edge(phage_1, phage_2, weight=row.jaccard_index)

        # dummy = g.es[0].attributes()

        g.vs["label"] = g.vs["name"]
        g.vs["color"] = "blue"

        # g.es["label"] = g.es["weight"]

        figure_name = "{}{}{}{}ph_plots.graph_for_cluster.{}.{}.pdf" \
            .format(self.mydir, self.dir_sep, "mcl_75.I25.plots", self.dir_sep, phc, self.extension2)

        # layout = g.layout_fruchterman_reingold()
        layout = g.layout_fruchterman_reingold(weights=g.es["weight"])
        ig.plot(g, figure_name, layout=layout)

phc_eval = PhageClusterEvaluation(mydir, extension, extension2)

phc_eval.read_files()

phc_eval.determine_origin()

phc_eval.make_jaccard_distribution()

phc_eval.make_cluster_graphs()
