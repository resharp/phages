import pandas as pd

class IctvStats :
    filename = ""
    ictv_data = None
    prok_data = None

    def __init__(self, filename):
        self.filename = filename
        self.ictv_data = pd.read_csv(filename, delimiter=";")

        # add derived fields
        self.ictv_data['count'] = 1

        species = self.ictv_data.Species
        hosts = []
        for row in species:
            hosts.append(row.split(" ")[0])
        self.ictv_data['host'] = hosts

    def nr_of_species(self):

        col = 'Species'
        species_names = self.ictv_data[col]
        unique_species_names = list(set(species_names))

        return len(unique_species_names)

    def nr_of_hosts(self):

        col = 'Species'
        species_names = self.ictv_data[col]
        hosts = [s.split(" ")[0] for s in species_names]

        unique_hosts = list(set(hosts))

        return len(unique_hosts)

    def nr_of_caudovirales_species(self):
        # cols = ictv_data.columns

        ## you can also use where but that changes just the complete row to NaN
        cau_data = self.ictv_data[self.ictv_data['Order'] == 'Caudovirales']

        return len(cau_data)

    def nr_of_non_blank_order(self):

        #comparing to itself is false for anything that is nan, which seems to be in it for the "empty" Order columns
        blank_order = self.ictv_data[self.ictv_data['Order'] != self.ictv_data['Order']]

        nr_blank_order = len(blank_order)
        nr_species = self.nr_of_species()
        nr_non_blank_order = nr_species - nr_blank_order

        ratio = float(nr_non_blank_order / nr_species)

        return ratio

    def nr_of_families(self):
        families = self.ictv_data['Family']

        unique_families = list(set(families))

        return len(unique_families)

    def top_five_hosts(self, order):

        # here me make a copy of the data
        cau_data = self.ictv_data[self.ictv_data['Order'] == order]

        top_counts = self.get_top_counts(cau_data, 'host')

        # cau_data[['host', 'count']].groupby('host').sum().sort_values(by='count', ascending=False)

        top_hosts = []
        for col in top_counts.axes[0][0:5]:
             top_hosts.append(col)

        return top_hosts

    # this might be a very complicated way, but it works
    def get_top_counts(self, df, col):
            return df[[col, 'count']].groupby(col).sum().sort_values(by='count', ascending=False)

    def top_hosts(self, order):
        return self.top_five_hosts(order)

    def print_top_host_counts(self):
        print(self.get_top_counts(self.ictv_data, 'host'))

    def print_top_family_counts(self):
        print(self.get_top_counts(self.ictv_data, 'Family'))

    def get_phages(self):
        #TODO inner join with other table prokaryotes.csv

        self.prok_data= pd.read_csv("prokaryotes.csv", delimiter=";")

        join_data = self.ictv_data.merge(self.prok_data, left_on=self.ictv_data.host, right_on=self.prok_data.host,how='inner')

        print(len(join_data))

        return join_data
