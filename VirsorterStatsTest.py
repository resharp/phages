import unittest

from MyTestCase import MyTestCase
from VirsorterStats import VirsorterStats
import configparser
import random
import logging

class VirsorterStatsTest(MyTestCase):

    config = configparser.ConfigParser()
    config.read('VirsorterStats.ini')

    mydir = config['output']['Directory']

    vs_stats = VirsorterStats(mydir)

    logging.basicConfig(filename='VirsorterStatsTest.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
    logging.info('Start VirsorterStatsTest')

    def test_vs_stats_non_empty(self):
        self.assertIsNotNone(self.vs_stats)

    def test_phage_data_non_empty(self):
        self.assertIsNotNone(self.vs_stats.phage_data)

    def test_phage_data_columns_contains_essential_fields(self):

        self.log_start()

        cols = self.vs_stats.phage_data.columns.values
        assert 'Category' in cols
        assert 'Contig_id' in cols
        assert 'phage_gene_start' in cols
        assert 'phage_gene_end' in cols


    def test_counts_of_category_1_is_larger_than_0(self):

        self.log_start()
        assert self.vs_stats.category_counts()['1'] > 0

    def test_nr_of_phages_larger_than_10(self):

        self.log_start()

        self.assertGreater(self.vs_stats.nr_of_phages(), 10)

    def test_nr_of_phages_equals_sum_category_counts(self):
        self.log_start()
        keys = ['1', '2', '3']
        total = 0
        for key in keys:
            total += self.vs_stats.category_counts()[key]
        self.assertEqual(total, self.vs_stats.nr_of_phages())

    def test_nr_of_complete_contigs_larger_than_1(self):
        self.log_start()
        self.assertGreater(self.vs_stats.nr_of_complete_contigs(), 1)

    def test_affi_contains_non_empty_genes(self):
        self.log_start()
        sub_list_empty_genes = self.vs_stats.all_affi_data[self.vs_stats.all_affi_data.gene_name == ""]

        self.assertEqual(len(sub_list_empty_genes), 0)

        sub_list_dash_genes = self.vs_stats.all_affi_data[self.vs_stats.all_affi_data.gene_name == "-"]
        self.assertEqual(len(sub_list_dash_genes), 0)

        # print("Length of filtered affi contig file:")
        # print(len(self.vs_stats.all_affi_data))

    #this method probably is going to confuse?
    #because we are not looking at the shared genes here
    #TODO: remove this test
    def test_top_gene_is_individual_gene(self):
        self.log_start()

        best_gene = self.vs_stats.all_gene_counts().index[0]

        logging.info("The most abundant gene in all contigs (not just phages)")
        logging.info(best_gene)

        self.assertEqual(best_gene.find("gi_"), 0)

    def test_convert_gene_id_to_contig_id(self):
        self.log_start()

        gene_id = "VIRSorter_PQXB01000001___PQXB01000001_1___[ANME-1_cluster_archaeon_strain_G37ANME1___2056316_3]-gene_8"
        contig_id = "VIRSorter_PQXB01000001___PQXB01000001_1___[ANME-1_cluster_archaeon_strain_G37ANME1___2056316_3]"

        self.assertEqual(self.vs_stats.gene_to_contig(gene_id), contig_id)

    def test_convert_gene_id_to_nr(self):
        self.log_start()

        gene_id = "VIRSorter_PQXB01000001___PQXB01000001_1___[ANME-1_cluster_archaeon_strain_G37ANME1___2056316_3]-gene_8"
        gene_nr = "8"
        self.assertEqual(self.vs_stats.gene_to_nr(gene_id), gene_nr)

    def test_phage_affi_data_columns_contains_essential_fields(self):
        self.log_start()

        cols = self.vs_stats.phage_affi_data.columns.values

        assert 'key_0' not in cols
        assert 'Contig_id' in cols
        assert 'gene_nr' in cols

    def test_phage_affi_contains_contigs_with_predicted_phages(self):

        self.log_start()

        phage_contigs = list(self.vs_stats.phage_data.Contig_id)

        list_of_contigs_in_affi = list(self.vs_stats.phage_affi_data.Contig_id.head(100))

        for i in range(0,5):
            random_affi_contig = random.choice(list_of_contigs_in_affi)
            # print(random_affi_contig)
            assert random_affi_contig in phage_contigs

    def test_most_abundant_shared_genes_contain_phage_clusters(self):

        self.log_start()

        #https://stackoverflow.com/questions/35268817/unique-combinations-of-values-in-selected-columns-in-pandas-data-frame-and-count
        combi_counts = self.vs_stats.contig_gene_counts()

        assert len(combi_counts) > 0

    def test_phage_affi_only_contains_genes_that_are_part_of_the_fragment(self):

        self.log_start()

        #example:
        #VIRSorter_MWYY01000037___MWYY01000037_1___[Acidobacteria_bacterium_28-9_strain_28-9___1962176_3]-gene_2-gene_16
        #gene_1 and gene_17 should not be included

        df = self.vs_stats.phage_affi_data
        unfiltered_phage_affi = df

        len_unfiltered = len(unfiltered_phage_affi)

        filtered_phage_affi = df[
            (df.gene_nr >= df.phage_gene_start)
            & (df.gene_nr <= df.phage_gene_end)
        ]

        len_filtered = len(filtered_phage_affi)

        self.assertEqual(len_unfiltered, len_filtered)

    #TODO look at the top Phage clusters and genes of the shared files
    def test_mostly_shared_genes_do_not_contain_genes_used_only_once(self):

        self.log_start()

        assert 1 == 1

    def test_mostly_shared_genes_contain_simple_annotation(self):
        self.log_start()

        assert 1 == 1

if __name__ == '__main__':
    unittest.main()
