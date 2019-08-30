import pandas as pd
import os
import logging

#TODO: set clustering result you want to evaluate
extension = "mcl_75.I20"

#TODO: set directory of predicted proteins and mcl clustering results
if os.name == "nt":
    eval_dir = "D:\\17 Dutihl Lab\\_tools\mcl"
else:
    eval_dir = "/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000"

class MclClusterEvaluation:

    extension = ""
    eval_dir = ""

    ip_pc_table_name = ""    # specific mcl result (clustering of IPs into PCs)
    pc_table_name = ""
    out_mcl_eval_name = ""

    ip_pc_df = None
    pc_df = None

    hit_df = None        # temp df containing HMM hits for one PC against all IPs (not just from one PC)
    ip_df = None         # temp df containing IPs for one PC

    mcl_eval_out_file = None                    # output file

    def __init__(self, eval_dir, extension):

        logging.basicConfig(filename='MclClusterEvaluation.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.DEBUG)

        self.eval_dir = eval_dir
        self.extension = extension

        if os.name == 'nt':
            dir_sep = "\\"
        else:
            dir_sep = "/"

        self.ip_pc_table_name = self.eval_dir + dir_sep + "ip_pc_table." + self.extension
        self.pc_table_name = self.eval_dir + dir_sep + "pc_table." + self.extension
        self.out_mcl_eval_name = self.eval_dir + dir_sep + "mcl_eval." + self.extension + ".txt"

    def print_file_names(self):

        print("Start evaluating mcl clusters for " + self.extension)
        print(self.ip_pc_table_name)
        print(self.pc_table_name)

    def read_files(self):

        logging.debug("start reading tables")
        self.read_ip_pc_table()
        self.read_pc_table()
        logging.debug("finished reading tables")

    def read_ip_pc_table(self):
        self.ip_pc_df = pd.read_csv(self.ip_pc_table_name, delimiter=" ",
                                    header=None
                                    , usecols=[0, 1]
                                    , names=['ip_id', 'pc_id'])

    def read_pc_table(self):
        self.pc_df = pd.read_csv(self.pc_table_name, delimiter=" ",
                                 header=None
                                 , usecols=[0]
                                 , names=['pc_id'])


    def eval_clusters(self):

        self.mcl_eval_out_file = open(self.out_mcl_eval_name, "w+")

        self.mcl_eval_out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    "pc_id", "max_f1_score", "tot_cl_ips", "tot_hits", "nr_tp", "nr_fp"
                                    , "high_bitscore"
                                    , "nc_bitscore", "nc_loc", "nc_ip"
                                    , "tc_bitscore", "tc_loc", "tc_ip"
                                    , "width_tc_nc", "width_tc_nc_score"))

        #TODO: use head if you want to process only a subset of top clusters
        #self.pc_df = self.pc_df.head(5)

        for index, row in self.pc_df.iterrows():
            pc_id = row.pc_id
            logging.debug("processing " + pc_id)

            self.eval_cluster(pc_id)

        self.mcl_eval_out_file.close()

    # added scores for recall, precision, and f1_score and calculate max f1_score for each PC
    # recall = sum_tps/tot_tps
    # precision = sum_tps/sum_hits
    def eval_cluster(self, pc_id):

        pc_dir = self.extension

        if os.name == 'nt':
            dir_sep = "\\"
        else:
            dir_sep = "/"

        hmm_result_table_name = self.eval_dir + dir_sep + pc_dir + dir_sep + pc_id + "_mafft_hmm_results_all_table.txt"
        one_pc_table_name = self.eval_dir + dir_sep + pc_dir + dir_sep + pc_id + ".txt"

        if not os.path.isfile(hmm_result_table_name):
            log_line = hmm_result_table_name + " does not exist."
            logging.warning(log_line)
            return

        log_line = "Evaluating " + hmm_result_table_name
        logging.debug(log_line)
        print(log_line)

        self.hit_df = pd.read_csv(hmm_result_table_name, comment="#"
                                  , delim_whitespace=True
                                  , header=None
                                  , usecols=[0,4,5]
                                  , names=["ip_hit", "evalue", "bitscore"] )

        self.ip_df = pd.read_csv(one_pc_table_name
                                ,header=None
                                ,usecols=[0]
                                ,names=["ip_cluster"])

        # print("---------------------------")
        # print("Number of IPs in cluster " + pc_id + ": " + str(len(self.ip_df)))

        merge_df = self.hit_df.merge(self.ip_df
                                     , left_on=self.hit_df.ip_hit
                                     , right_on=self.ip_df.ip_cluster
                                     , how='left')


        false_positives = merge_df[merge_df.isnull().ip_cluster]

        true_positives = merge_df[merge_df.ip_cluster.notnull()]

        total_tps = len(true_positives)
        sum_tps = 0

        f1_scores = []
        for index, row in merge_df.iterrows():

            sum_hits = index + 1
            if str(row.ip_cluster) != "nan":
                sum_tps += 1

            recall = sum_tps / total_tps
            precision = sum_tps / sum_hits

            if recall + precision > 0:
                f1_score = 2 * recall * precision / (recall + precision)
            else:
                f1_score = 0

            f1_scores += [f1_score]

        max_f1_score = max(f1_scores)

        tot_cl_ips = str(len(self.ip_df))
        tot_hits = str(len(merge_df))
        nr_fp = str(len(false_positives))
        nr_tp = str(total_tps)

        high_bitscore = str(merge_df.iloc[0].bitscore)

        #noise cutoff, only if false positives exist
        if len(false_positives) > 0:
                index_nc_1 = false_positives.index[0]

                nc_bitscore = str(merge_df.iloc[index_nc_1].bitscore)
                nc_loc = str(index_nc_1 + 1)
                nc_ip = str(false_positives.iloc[0].ip_hit)
        else:
                nc_bitscore = "-"
                nc_loc = "-"
                nc_ip = "-"

        #trusted cutoff, only if true positives exist
        if len(true_positives) > 0:

                index_tc_1 = true_positives.index[-1]

                tc_bitscore = str(merge_df.iloc[index_tc_1].bitscore)
                tc_loc = str(index_tc_1 + 1)
                tc_ip = true_positives.iloc[-1].ip_hit
        else:
                tc_bitscore = "-"
                tc_loc = "-"
                tc_ip = "-"

        #width between nc and tc
        if len(true_positives) > 0 and len(false_positives) > 0:
                width_tc_nc = str(index_tc_1 - index_nc_1)
                width_tc_nc_score = str(merge_df.iloc[index_nc_1].bitscore - merge_df.iloc[index_tc_1].bitscore)
        else:
                width_tc_nc = "-"
                width_tc_nc_score = "-"

        self.mcl_eval_out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    pc_id, max_f1_score, tot_cl_ips, tot_hits, nr_tp, nr_fp
                                    , high_bitscore
                                    , nc_bitscore, nc_loc, nc_ip
                                    , tc_bitscore, tc_loc, tc_ip
                                    , width_tc_nc, width_tc_nc_score ))

cluster_eval = MclClusterEvaluation(eval_dir, extension)

cluster_eval.print_file_names()

cluster_eval.read_files()

cluster_eval.eval_clusters()
