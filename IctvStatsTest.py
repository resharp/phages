import unittest

from IctvStats import IctvStats


#TODO show how to run test suites
class IctvStatsTest(unittest.TestCase) :

    mydir = "D:\\17 Dutihl Lab\\source\\phages_pycharm\\"

    ictv_stats = IctvStats(mydir)

    def test_filename(self):
        self.assertTrue(self.ictv_stats.filename != "")

    def test_ictv_data_non_empty(self):
        self.assertIsNotNone(self.ictv_stats.ictv_data)

    def test_nr_species_larger_than_100(self):
        self.assertGreater(self.ictv_stats.nr_of_species(), 100)

    #test number of hosts is larger than 100
    def test_nr_hosts_larger_than_100(self):
        self.assertGreater(self.ictv_stats.nr_of_hosts(), 100)

    #test number of Caudovirales is larger than zero
    def test_nr_of_caudovirales_is_1320(self):
        self.assertEqual(self.ictv_stats.nr_of_caudovirales_species(), 1320)

    #test percentage of non-blank Order is larger than 50%
    def test_nr_of_non_blank_order_is_larger_than_50_percent(self):
        nr_non_blank_order = self.ictv_stats.nr_of_non_blank_order()
        self.assertGreater(nr_non_blank_order, 0.50)

    def test_nr_of_families_larger_than_100(self):
        self.assertGreater(self.ictv_stats.nr_of_families(), 100)

    def test_top_five_hosts_for_Caudovirales_contain_Escherichia(self):

        species = "Escherichia";

        assert species in self.ictv_stats.top_five_hosts('Caudovirales'), format("{0} should be in top five hosts", species)

    def test_top_hosts_default_length_five(self):
        self.assertEqual(len(self.ictv_stats.top_hosts('Caudovirales')), 5)

    # def test_top_five_hosts(self):
    #     print(self.ictv_stats.top_five_hosts('Caudovirales'))

    # def test_print_top_family_counts(self):
    #     self.ictv_stats.print_top_family_counts()

    def test_nr_of_phages_smaller_than_total(self):

        nr_of_phages = len(self.ictv_stats.get_phages())
        self.assertLess(nr_of_phages, 5000)

## initialization of test fixture
## all the initialization goes in here to keep the tests super focused
if __name__ == '__main__':

    unittest.main()