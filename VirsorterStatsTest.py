import unittest

from VirsorterStats import VirsorterStats
import configparser
import random

class VirsorterStatsTest(unittest.TestCase):

    config = configparser.ConfigParser()
    config.read('VirsorterStats.ini')

    mydir = config['output']['Directory']

    vs_stats = VirsorterStats(mydir)

    def test_vs_stats_non_empty(self):
        self.assertIsNotNone(self.vs_stats)

    def test_phage_data_non_empty(self):
        self.assertIsNotNone(self.vs_stats.phage_data)

    def test_phage_data_columns_contains_essential_fields(self):
        cols = self.vs_stats.phage_data.columns.values
        assert 'Category' in cols
        assert 'Contig_id' in cols

    def test_counts_of_category_1_is_2(self):
        assert self.vs_stats.category_counts()['1'] == 2

    def test_nr_of_phages_larger_than_10(self):
        self.assertGreater(self.vs_stats.nr_of_phages(), 10)

    def test_nr_of_phages_equals_sum_category_counts(self):
        keys = ['1', '2', '3']
        total = 0
        for key in keys:
            total += self.vs_stats.category_counts()[key]
        self.assertEqual(total, self.vs_stats.nr_of_phages())

    def test_nr_of_complete_contigs_larger_than_1(self):
        self.assertGreater(self.vs_stats.nr_of_complete_contigs(), 1)

    def test_affi_contains_non_empty_genes(self):

        sub_list_empty_genes = self.vs_stats.all_affi_data[self.vs_stats.all_affi_data.gene_name == ""]

        self.assertEqual(len(sub_list_empty_genes), 0)

        sub_list_dash_genes = self.vs_stats.all_affi_data[self.vs_stats.all_affi_data.gene_name == "-"]
        self.assertEqual(len(sub_list_dash_genes), 0)

        # print("Length of filtered affi contig file:")
        # print(len(self.vs_stats.all_affi_data))

    #this method probably is going to confuse
    #because we are not looking at the shared genes (yet)
    def test_top_gene_is_individual_gene(self):
        best_gene = self.vs_stats.all_gene_counts().index[0]
        print("The most abundant gene in all contigs (not just phages)")
        print(best_gene)
        self.assertEqual(best_gene.find("gi_"), 0)

    # def test_convert_gene_id_to_contig_id(self):
    #
    #     gene_id = "VIRSorter_PQXB01000001___PQXB01000001_1___[ANME-1_cluster_archaeon_strain_G37ANME1___2056316_3]-gene_8"
    #     contig_id = "VIRSorter_PQXB01000001___PQXB01000001_1___[ANME-1_cluster_archaeon_strain_G37ANME1___2056316_3]"
    #
    #     self.assertEqual(self.vs_stats.gene_to_contig(gene_id), contig_id)

    def test_phage_affi_data_columns_contains_essential_fields(self):
        cols = self.vs_stats.phage_affi_data.columns.values

        assert 'key_0' not in cols
        assert 'Contig_id' in cols

    def test_phage_affi_contains_contigs_with_predicted_phages(self):

        phage_contigs = list(self.vs_stats.phage_data.Contig_id)

        list_of_contigs_in_affi = list(self.vs_stats.phage_affi_data.Contig_id)

        for i in range(0,5):
            random_affi_contig = random.choice(list_of_contigs_in_affi)
            # print(random_affi_contig)
            assert random_affi_contig in phage_contigs

# also join the simple annotations from the affi file
# look at the top Phage clusters and genes of the shared files

if __name__ == '__main__':
    unittest.main()
