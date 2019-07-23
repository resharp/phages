import pandas as pd
import os

# self.phage_affi_data[self.phage_affi_data.gene_name.str.find("Phage_cluster") == 0].gene_name.value_counts()

class VirsorterStats:
    phage_signal_filename = ""
    phage_affi_filename = ""
    phage_data = None
    phage_affi_data = None

    def __init__(self, mydir):

        self.phage_signal_filename = mydir + "VIRSorter_global-phage-signal.csv"
        self.phage_affi_filename = mydir + "Metric_files\\VIRSorter_affi-contigs.tab"

        if not os.path.isfile(self.phage_affi_filename):
            #switch to linux style
            self.phage_affi_filename = mydir + "Metric_files/VIRSorter_affi-contigs.tab"

        self.read_phage_data(self.phage_signal_filename)
        self.read_affi_data(self.phage_affi_filename)

    def read_phage_data(self, phage_signal_filename):

        f = open(phage_signal_filename, 'r')
        lines = f.readlines()
        header = lines[1].replace("## ", "")
        lines_no_comment = [n for n in lines if not n.startswith('##')]
        f.close()

        filename2 = phage_signal_filename + "2"
        f2 = open(filename2, 'w')
        f2.write(header)
        f2.writelines(lines_no_comment)
        f2.close()

        self.phage_data= pd.read_csv(filename2, delimiter=",")

    def read_affi_data(self, phage_affi_filename):

        f = open(phage_affi_filename, 'r')
        lines = f.readlines()
        lines_no_comment = [n for n in lines if not n.startswith('>')]
        f.close()

        filename2 = phage_affi_filename + "2"
        f2 = open(filename2, 'w')
        f2.write("gene_id|nr1|nr2|nr3|nr4|gene_name|nr_6|nr_7|nr_8|short_annotation|nr_10|nr_11")
        f2.writelines(lines_no_comment)
        f2.close()

        self.phage_affi_data = pd.read_csv(filename2, delimiter="|")

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





