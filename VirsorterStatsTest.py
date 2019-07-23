import unittest

from VirsorterStats import VirsorterStats


class VirsorterStatsTest(unittest.TestCase):

    mydir = "D:\\17 Dutihl Lab\\virsorter_output\\vs_all_100\\"
    vs_stats = VirsorterStats(mydir)

    def test_vs_stats_non_empty(self):
        self.assertIsNotNone(self.vs_stats)

    def test_phage_data_non_empty(self):
        self.assertIsNotNone(self.vs_stats.phage_data)

    def test_phage_data_columns_contains_category(self):
        cols = self.vs_stats.phage_data.columns.values
        assert 'Category' in cols

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

# test that the affi only contains non-empty genes
# we are only interested in the annotations
# also join the simple annotations from the affi file
# look at the top Phage clusters and genes of the shared files

# extract contigs

if __name__ == '__main__':
    unittest.main()
