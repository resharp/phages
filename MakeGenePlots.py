import argparse
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns;sns.set()
import sys

# MakeGenePlots
# create multiple heatmaps
#   1st with the log pN/pS values
#   2nd with the missing genes (just 0/1 for clarity)
#   3rd distribution plots [integration over structural genes or other categories]
#
#   all sample_gene measures are in small separate files
#   we aggregate them to gene level (aggregation over all samples)
#   therefore we load all files and merge them into one dataframe
#
# masking missing values in heatmaps:
# https://github.com/mwaskom/seaborn/issues/375


class MakeGenePlots:

    sample_dir = ""
    ref_dir = ""
    dir_sep = ""
    plot_dir = ""
    ref = ""  # temporary, it may be that MakeGenePlots will go beyond the level of ref

    gene_sample_df = None
    filtered_gene_sample_df = None
    gene_df = None
    sample_df = None
    gene_anno_df = None

    bin_sample_df = None
    bin_df = None

    sample_breadth_df = None

    threshold_depth = 0
    threshold_breadth = 0

    all_scores_df = None

    def __init__(self, sample_dir, ref_dir, ref, threshold_depth, threshold_breadth):

        self.sample_dir = sample_dir
        self.ref_dir = ref_dir
        self.ref = ref
        self.threshold_depth = threshold_depth
        self.threshold_breadth = threshold_breadth

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=self.sample_dir + self.dir_sep + 'MakeGenePlots.log', filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.INFO)

    def read_files_for_ref(self, ref):

        logging.info("start reading tables")

        self.gene_sample_df = self.read_and_concat_measures("_gene_measures.txt", ref)

        self.gene_anno_df = self.read_gene_annotation(ref)

        self.bin_sample_df = self.read_and_concat_measures("_bin_measures.txt", ref)

        self.sample_breadth_df = self.read_and_concat_measures("_sample_measures.txt", ref)

        logging.info("finished reading tables")

    def read_and_concat_measures(self, file_name_ext, ref=None):

        if ref:
            subfolders = [f.path for f in os.scandir(self.sample_dir) if f.is_dir() and ref in f.name]
        else:
            subfolders = [f.path for f in os.scandir(self.sample_dir) if f.is_dir() and "GenePlots" not in f.name]

        samples = []

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            sample = sample.split("_")[0]
            sample_name = subfolder + self.dir_sep + sample + file_name_ext

            if os.path.isfile(sample_name) and os.path.getsize(sample_name) > 0:
                logging.info("processing {}".format(sample_name))

                gene_sample_df = pd.read_csv(sample_name
                                         ,  sep='\t'
                                         )
                samples.append(gene_sample_df)

        return pd.concat(samples)

    def read_gene_annotation(self, ref):

        gene_anno_file_name = self.ref_dir + self.dir_sep + "{ref}_gene_list.txt".format(ref=ref)

        # to do also read gene_fam_x (to be produced by AnnotateCrassGenomes)
        anno_df = pd.read_csv(gene_anno_file_name
                              ,   sep='\t'
                              ,   header=None
                              ,   usecols=[0, 1, 2, 3]
                              ,   names=["Protein", "gene_fam", "region", "Annotation"]
                              )
        return anno_df

    def prepare_data_for_plots(self):

        # here we can prepare the data further for displaying purposes
        # it could also have been done in CalcDiversiMeasures

        # we want to see the missing genes, and missing genes are if there are no mappings or hardly any
        self.gene_sample_df['missing_gene'] = 0
        self.gene_sample_df.loc[self.gene_sample_df.AAcoverage_perc < 0.02 , 'missing_gene'] = 1
        self.gene_sample_df.loc[np.isnan(self.gene_sample_df.AAcoverage_perc), 'missing_gene'] = 1

        self.gene_sample_df['double_coverage'] = 0
        self.gene_sample_df.loc[self.gene_sample_df.AAcoverage_perc > 2 , 'double_coverage'] = 1

        # data_debug = self.gene_df[['sample','Protein', 'AAcoverage_perc', 'missing_gene', 'double_coverage']]

        # make binary plot for positive selection or conservation
        self.gene_sample_df['positive_selection'] = 0
        self.gene_sample_df.loc[self.gene_sample_df["log10_pN/pS"] > -0.3, 'positive_selection'] = 1
        self.gene_sample_df.loc[self.gene_sample_df["log10_pN/pS"] < -0.7, 'positive_selection'] = -1

    def create_plot_dir(self, ref):

        self.plot_dir = self.sample_dir + self.dir_sep + ref + "_" + "GenePlots"
        os.makedirs(self.plot_dir, exist_ok=True)

    def create_histograms(self):

        data = self.filtered_gene_sample_df

        # we only need the counts per gene
        counts_df = data.groupby("Protein").agg(
            {
                'log10_pN/pS': ["count"]
            }).reset_index()
        counts_df.columns = ["_".join(x) for x in counts_df.columns.ravel()]
        counts_df.rename(columns={'Protein_': 'Protein'}, inplace=True)

        min = counts_df['log10_pN/pS_count'].min()
        max = counts_df['log10_pN/pS_count'].max()

        sns.distplot(counts_df[['log10_pN/pS_count']], bins=(max-min+1), kde=False)

        plt.title("{ref}: distribution of sample presence for genes ({breadth}/{depth}x)".
                  format(ref=self.ref, depth=self.threshold_depth, breadth=self.threshold_breadth))
        figure_name = "{}{}gene_plots.gene_sample_counts.{breadth}.{depth}x.pdf".\
            format(self.plot_dir, self.dir_sep, depth=self.threshold_depth, breadth=self.threshold_breadth)

        plt.savefig(figure_name)

    def create_unfiltered_heat_maps(self):

        # use unfiltered data for quality measures (and do not add suffix)
        data = self.gene_sample_df

        self.create_heat_map(data, "AAcoverage_cv", "Internal coefficient of variation per gene", add_suffix=False)
        self.create_heat_map(data, "AAcoverage_perc", "Coverage % (compared to whole genome)", add_suffix=False)

        self.create_heat_map(data, "breadth_1x", "1x horizontal coverage percentages", add_suffix=False)
        self.create_heat_map(data, "breadth_10x", "10x horizontal coverage percentages", add_suffix=False)
        self.create_heat_map(data, "breadth_50x", "50x horizontal coverage percentages", add_suffix=False)

        self.create_heat_map(data, "missing_gene"
                             , "Missing genes (based on coverage < 2% rest genome)", add_suffix= False)

        self.create_heat_map(data, "double_coverage",
                             "Genes that have > twice the amount of coverage compared to genome", add_suffix=False)

    def create_filtered_heat_maps(self):

        # also filter out complete samples before showing pN/pS and other measures
        # we do not care for samples with a 1x horizontal coverage below 5%
        filtered_data = self.filtered_gene_sample_df

        # make a heat map of the log10_pN/pS based on multiple samples
        self.create_heat_map(filtered_data, "log10_pN/pS", "Log 10 of pN/pS (red=diverging, blue=conserved)")
        self.create_heat_map(filtered_data, "positive_selection", "log10_pN/pS either > -0.3 or < - 0.7")

        self.create_heat_map(filtered_data, "SndAAcnt_perc_polymorphism_mean", "Within sample AA variation in genes")

        self.create_heat_map(filtered_data, "entropy_mean", "Mean codon entropy (base 10)")

    def filter_gene_sample_on_sample_and_gene_coverage(self, filter_genes=False):

        data = self.gene_sample_df

        filtered_sample = self.sample_breadth_df[self.sample_breadth_df.breadth_1x.ge(0.05)]
        filtered_sample = filtered_sample.drop(["breadth_1x", "breadth_10x", "breadth_50x", "breadth_100x"], axis=1)

        data = data.merge(filtered_sample,
                          left_on=data["sample"],
                          right_on=filtered_sample["sample"],
                          how="inner").drop(["key_0"], axis=1).drop(["sample_y"], axis=1)

        data.rename(columns={'sample_x': 'sample'}, inplace=True)

        if filter_genes:
            breadth_field = "breadth_{depth}x".format(depth=self.threshold_depth)
            data = data[data[breadth_field] > self.threshold_breadth]

        # add annotation, e.g. region
        merge_df = data.merge(self.gene_anno_df
                              , left_on=data.Protein
                              , right_on=self.gene_anno_df.Protein
                              , how='left').drop(["Protein_x", "Protein_y", "ref_y"], axis=1)
        merge_df.rename(columns={'key_0': 'Protein', "ref_x": "ref"}, inplace=True)

        # to do this assertion does not seem to hold: how come?
        # assert(len(data) == len(merge_df))

        self.filtered_gene_sample_df = merge_df

        # this was the old way of filtering (when we always had enough depth in the pilot):
        # data = data[(data.AAcoverage_perc < 1.5)]
        # data = data[(data.AAcoverage_perc > 0.2)]

    # https://seaborn.pydata.org/generated/seaborn.heatmap.html
    def create_heat_map(self, data, measure, title, add_suffix=True):

        data = data[['Protein', 'sample', measure]]
        data = data.sort_values("Protein")
        data = data.set_index(["sample", "Protein"])

        # nr_of_samples = 3
        # nr_of_genes = 10
        # data = data.tail(nr_of_samples * nr_of_genes)

        # convert from multi-index to cross-product table
        data = data.unstack()

        # data = data.transpose() # if you would like to switch columns and rows

        # rename columns, unstack
        data.columns = ["_".join(x) for x in data.columns.ravel()]

        # data.columns = [x.split("_")[-1] for x in data.columns]
        data.columns = ["_".join(x.split("_")[3:]) for x in data.columns]

        plt.clf()  # clear the current figure (always do this before any new plot)
        ax = sns.heatmap(data, cmap="seismic", annot=False)

        suffix_title = ""
        if add_suffix:
            suffix_title = " ({breadth}/{depth}x)".format(depth=self.threshold_depth, breadth=self.threshold_breadth)

        plt.title("{ref}: {title}{suffix}".format(ref=self.ref, title=title, suffix=suffix_title))

        suffix = ""
        if add_suffix:
            suffix = ".{breadth}.{depth}x".format(depth=self.threshold_depth, breadth=self.threshold_breadth)

        figure_name = "{}{}gene_plots.heat_map.{measure}{suffix}.pdf".format(
            self.plot_dir, self.dir_sep, measure=measure.replace("/", "_"),
            suffix=suffix
        )
        plt.savefig(figure_name)

    # we want to see what genes have the highest and lowest pN/pS scores
    # based on self.gene_sample_df
    # that can be further aggregated into self.gene_df
    # and then merged with self.gene_annotation.df for annotation
    def score_genes_and_output_ordered_genes_to_file(self):

        self.gene_df = self.filtered_gene_sample_df.groupby("Protein").agg(
            {
                'log10_pN/pS': ["mean", "count", "std"]
            ,   'entropy_mean': ["mean", "count"]
            }).reset_index()

        self.gene_df.columns = ["_".join(x) for x in self.gene_df.columns.ravel()]
        self.gene_df.rename(columns={'Protein_': 'Protein'}, inplace=True)

        self.gene_df = self.gene_df.sort_values(by='log10_pN/pS_mean', ascending=False)

        merge_df = self.gene_df.merge(self.gene_anno_df
                                      , left_on=self.gene_df.Protein
                                      , right_on=self.gene_anno_df.Protein
                                      , how='inner')
        merge_df.rename(columns={'key_0': 'Protein'}, inplace=True)

        # to do: also add gene_fam as a key for combining Protein later on
        merge_df = merge_df[['Protein', 'gene_fam',
                             'log10_pN/pS_mean', 'log10_pN/pS_std', 'log10_pN/pS_count', 'Annotation']]

        filename = self.plot_dir + self.dir_sep + "crassphage_pN_pS_values.{breadth}.{depth}x.txt".format(
            depth=self.threshold_depth, breadth=self.threshold_breadth
        )
        merge_df.to_csv(filename, index=False, sep='\t')

    def score_samples(self):

        self.sample_df = self.filtered_gene_sample_df.groupby("sample").agg(
            {
                'log10_pN/pS':      ["mean", "count", "std"]
            ,   'entropy_mean':     ["mean", "count"]
            ,   'AAcoverage_mean':  ["mean"]
            }).reset_index()

        self.sample_df.columns = ["_".join(x) for x in self.sample_df.columns.ravel()]
        self.sample_df.rename(columns={'sample_': 'sample'}, inplace=True)

        self.sample_df = self.sample_df.sort_values(by='entropy_mean_mean', ascending=False)

    def create_box_plots(self, min_nr_samples=3):

        # box plots for log10_pN/pS

        # filter out the genes that do not have a minimal presence of min_nr_of_samples occurences in the samples
        filter_data = self.gene_df[self.gene_df['log10_pN/pS_count'] >= min_nr_samples]

        # head only works because self.gene_df is already ordered by log10_pN/pS_mean in score_genes_..()
        top10_data = filter_data.head(10)[['Protein', 'log10_pN/pS_mean']]
        self.create_box_plot(top10_data, "Protein", "log10_pN/pS",
                             "top 10 pos selection present in at least {} samples".format(min_nr_samples))

        bottom10_data = filter_data.tail(10)[['Protein', 'log10_pN/pS_mean']]
        self.create_box_plot(bottom10_data, "Protein", "log10_pN/pS",
                             "top 10 most conserved present in at least {} samples".format(min_nr_samples))

        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "Protein", "log10_pN/pS",
                             "top and bottom 10 present in at least {} samples".format(min_nr_samples))

        self.create_box_plot(self.gene_df, "Protein", "log10_pN/pS", "all genes")

        # box plots for ENTROPY
        filter_data = self.gene_df[self.gene_df['entropy_mean_count'] >= min_nr_samples]

        top10_data = filter_data.head(10)[['Protein', 'entropy_mean_mean']]
        self.create_box_plot(top10_data, "Protein", "entropy_mean",
                             "top 10 internal var present in at least {} samples".format(min_nr_samples))

        bottom10_data = filter_data.tail(10)[['Protein', 'entropy_mean_mean']]
        self.create_box_plot(bottom10_data, "Protein", "entropy_mean",
                             "bottom 10 internal var present in at least {} samples".format(min_nr_samples))

        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "Protein", "entropy_mean",
                             "top and bottom 10 internal var in at least {} samples".format(min_nr_samples))

        # new aggregation per sample (min_nr_samples should now be read as min_nr_genes)
        filter_data = self.sample_df[self.sample_df['log10_pN/pS_count'] >= min_nr_samples]

        top10_data = filter_data.head(10)
        bottom10_data = filter_data.tail(10)
        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "sample", "entropy_mean",
                             "top and bottom 10 entropy with at least {} genes".format(min_nr_samples))

        filter_data = filter_data.sort_values(by='AAcoverage_mean_mean', ascending=False)

        top10_data = filter_data.head(10)
        self.create_box_plot(top10_data, "sample", "AAcoverage_mean",
                             "top 10 coverage with at least {} genes".format(min_nr_samples))

        bottom10_data = filter_data.tail(10)
        self.create_box_plot(bottom10_data, "sample", "AAcoverage_mean",
                             "bottom 10 coverage with at least {} genes".format(min_nr_samples))

        combined_data = pd.DataFrame.append(top10_data, bottom10_data)
        self.create_box_plot(combined_data, "sample", "AAcoverage_mean",
                             "top and bottom 10 coverage with at least {} genes".format(min_nr_samples))

    def create_box_plot(self, filter_data, agg_field, measure, title):

        data = self.filtered_gene_sample_df

        if len(filter_data) == 0:
            logging.warning("no data left after filter applied: cannot create box plot "
                            "for {measure} aggregated on {agg_field}".format(measure=measure, agg_field=agg_field))
            return

        data = data.merge(filter_data
                             , left_on=data[agg_field]
                             , right_on=filter_data[agg_field]
                             , how='inner')

        data.rename(columns={"{}_x".format(agg_field): agg_field}, inplace=True)

        data = data[[agg_field, measure]]

        # add mean of measure per gene and then merge with the original dataset
        # in order to sort the genes on the mean of the measure
        grouped = data.groupby(agg_field).mean()

        data = data.merge(grouped
                          , left_on=data[agg_field]
                          , right_on = grouped.index
                          , how='inner').sort_values(["{}_y".format(measure), agg_field], ascending=False).\
            rename(columns={"{}_x".format(measure): measure })

        plt.clf()
        plt.title("{ref}: {title} ({breadth}/{depth}x)".
                  format(ref=self.ref, title=title, depth=self.threshold_depth, breadth=self.threshold_breadth))

        sns.set(style="ticks")

        # to do: We get a warning on the percentile calculations (implicit in box plot) for the infinite values
        # we should probably recalculate p_N/p_S with a pseudocount
        sns.boxplot(x=measure, y=agg_field, data=data,
                    whis="range", palette="vlag")

        sns.swarmplot(x=measure, y=agg_field, data=data,
                      size=2, color=".3", linewidth=0)

        figure_name = "{}{}gene_plots.{}.box_plot.{measure}.{title}.{breadth}.{depth}x.pdf".format(
            self.plot_dir, self.dir_sep, agg_field,
            measure=measure.replace("/", "_"), title=title.replace(" ", "_"),
            depth=self.threshold_depth, breadth=self.threshold_breadth
        )
        plt.savefig(figure_name)

    def create_region_plot(self):

        data = self.filtered_gene_sample_df
        data.loc[data.region == "assembly", "region"] = "assembly.rest"

        measure = "log10_pN/pS"
        agg_field = "region"
        title = "{measure} per region".format(measure=measure)

        plt.figure(figsize=(12, 5))
        plt.title("{ref}: {title} ({breadth}/{depth}x)".
                  format(ref=self.ref, title=title, depth=self.threshold_depth, breadth=self.threshold_breadth))

        sns.set(style="ticks")

        # to do: We get a warning on the percentile calculations (implicit in box plot) for the infinite values
        # we should probably recalculate p_N/p_S with a pseudocount
        sns.boxplot(x=measure, y=agg_field, data=data,
                    whis="range", color="white")

        sns.swarmplot(x=measure, y=agg_field, data=data,
                      size=2, color=".3", linewidth=0)

        figure_name = "{}{}gene_plots.{}.box_plot.{measure}.{title}.{breadth}.{depth}x.svg".format(
            self.plot_dir, self.dir_sep, agg_field,
            measure=measure.replace("/", "_"), title=title.replace(" ", "_").replace("/", "_"),
            depth=self.threshold_depth, breadth=self.threshold_breadth
        )

        plt.savefig(figure_name)

    def create_bin_line_plots_for_pn_ps(self):
        # make line plots based on
        # aggregation of bin_sample_df
        # to bin_df

        nr_lines = len(self.bin_sample_df)

        # to do: also change filter here (should also be changed in CalcDiversiMeasures)
        filtered_bin_data = \
            self.bin_sample_df[(self.bin_sample_df.AAcoverage_perc > 0.2) & (self.bin_sample_df.AAcoverage_cv < 0.2)]

        nr_filtered_lines = len(filtered_bin_data)

        # to do: do not take mean of log but instead log of mean
        data = filtered_bin_data.groupby(["Protein", "AAPosition"]).agg(
            {
                'log10_pN_pS_60': ["mean"]
            }).reset_index()
        data.columns = ["_".join(x) for x in data.columns.ravel()]
        data.rename(columns={'Protein_': 'Protein'}, inplace=True)
        data.rename(columns={'AAPosition_': 'AAPosition'}, inplace=True)

        data.rename(columns={'log10_pN_pS_60_mean': 'log10_pN_pS_60'}, inplace=True)

        data = data.join(data.groupby('Protein')['AAPosition'].max(), on='Protein', rsuffix='_max')

        self.line_plots_for_measure(data, "log10_pN_pS_60")

    def line_plots_for_measure(self, data, measure, ylim_bottom=None):

        genes = data.Protein.unique()

        # determine the number of genes and make a plot for every set of 6 genes
        nr_plots = int(np.round((len(genes)/6)))

        # we loop through the plots and subplots and then select the next gene instead of looping through the genes
        i = 0
        # we always make one extra plot (so some of the sub plots of the last plot may be empty)
        for j in range(0,nr_plots):

            fig, axs = plt.subplots(3,2)
            fig.tight_layout(pad=2)
            for row in axs:
                for ax in row:
                    #to prevent
                    if len(genes) > i:
                        gene = genes[i]
                        i = i + 1
                        gene_data = data[data.Protein == gene]

                        x_max = gene_data.AAPosition_max.max()

                        sns.lineplot(legend=None, data=gene_data, x="AAPosition", y=measure, ax=ax)
                        ax.set_ylim(bottom=ylim_bottom)
                        ax.set_xlim(left=1)
                        ax.set_xlim(right=x_max)
                        ax.set_title(gene)

            start = 1 + j*6
            end = 6 + j*6
            plot_name = self.plot_dir + self.dir_sep + "{}_{}_{}.png".format(measure, start, end)

            plt.savefig(plot_name)
            plt.close()

    def do_analysis(self, min_nr_samples, ref):

        self.read_files_for_ref(ref)

        self.create_plot_dir(ref)

        self.create_bin_line_plots_for_pn_ps()

        self.prepare_data_for_plots()

        self.create_unfiltered_heat_maps()

        self.filter_gene_sample_on_sample_and_gene_coverage(filter_genes=True)

        self.create_histograms()

        self.create_filtered_heat_maps()

        self.score_genes_and_output_ordered_genes_to_file()

        self.score_samples()

        # create box plots for top 10 and bottom 10 pN/pS values for genes with a minimum nr of samples for that gene
        # after applying the filter (data have been prepared in score_genes_..() and score_samples()
        self.create_box_plots(min_nr_samples)

        self.create_region_plot()

    #region family analysis
    def read_all_files_for_family(self):

        logging.info("start reading tables")

        self.gene_sample_df = self.read_and_concat_measures("_gene_measures.txt")

        # to do: read all annotations of all refs
        self.gene_anno_df = self.read_all_annotations()

        self.sample_breadth_df = self.read_and_concat_measures("_sample_measures.txt")

        logging.info("finished reading tables")

    def read_all_annotations(self):

        files = [f.path for f in os.scandir(self.ref_dir) if not f.is_dir() and "_gene_list.txt" in f.name]

        anno_list = []

        for file in files:
            anno_df = pd.read_csv(file
                                  , sep='\t'
                                  , header=None
                                  , usecols=[0, 1, 2, 3]
                                  , names=["Protein", "gene_fam", "region", "Annotation"]
                                  , skiprows=1
                                  )
            ref = file.split(self.dir_sep)[-1].replace("_gene_list.txt", "")
            anno_df["ref"] = ref

            anno_list.append(anno_df)

        return pd.concat(anno_list)

    def breadth_statistics(self):

        data = self.filtered_gene_sample_df

        refs = data.ref.unique()

        row_list = []

        breadth_field = "breadth_{depth}x".format(depth=self.threshold_depth)

        for ref in refs:
            data_ref = data[data.ref == ref]
            nr_data = len(data_ref)

            for threshold in range(1, 101):
                percentage = threshold/100
                nr_filtered = len(data_ref[data_ref[breadth_field] >= percentage])

                fraction = nr_filtered / nr_data

                row_list.append([percentage, ref, nr_filtered, fraction])

                pass

        df = pd.DataFrame(columns=("percentage", "ref", "count", "fraction"), data=row_list)

        self.make_cumulation_plot(df, refs, "fraction")

        self.make_cumulation_plot(df, refs, "count")

    def make_cumulation_plot(self, df, refs, measure):

        colors = sns.color_palette("husl", n_colors=len(refs))

        fig, ax = plt.subplots(figsize=(12, 5))

        ref_nr = 0

        if measure == "count":
            max_y = df[measure].max() + 1
        else:
            max_y = df[measure].max() + 0.01

        for ref in refs:

            df_ref = df[df.ref == ref]

            plt.xlim(0, 1)
            plt.ylim(0, max_y)

            sns.lineplot(x=df_ref.percentage,
                         y=df_ref[measure],
                         color=colors[ref_nr],
                         ax=ax)
            ref_nr = ref_nr + 1

        title = "{measure} of data points for genes for {depth}x breadth thresholds".format(
            depth=self.threshold_depth, measure=measure)
        plt.title(title)
        ax.legend(refs, facecolor='w')
        ax.set(xlabel='breadth percentage for gene'
               , ylabel='{measure} of gene data points used'.format(measure=measure))
        # plt.show()

        self.plot_dir = self.sample_dir + self.dir_sep + "FamilyPlots"
        figure_name = "{}{}family_plots.gene_usage.{measure}.{title}.{depth}x.svg".format(
            self.plot_dir, self.dir_sep,
            title=title.replace(" ", "_"),
            depth=self.threshold_depth, breadth=self.threshold_breadth,
            measure=measure
        )

        plt.savefig(figure_name)
        plt.clf()

    def rank_gene_families(self):

        data = self.filtered_gene_sample_df

        data = data[data.gene_fam.isnull() == False]

        # now first aggregate on gene_fam log10_pN/pS

        # to do: write this to text file
        df_fam = data.groupby("gene_fam").agg(
            {
                'log10_pN/pS': ["mean", "count"]
            }).reset_index()

        df_fam.columns = ["_".join(x) for x in df_fam.columns.ravel()]
        df_fam.rename(columns={'gene_fam_': 'gene_fam'}, inplace=True)

        df_fam = df_fam.sort_values(by='log10_pN/pS_mean', ascending=False)

        data_new = data.merge(df_fam,
                              left_on=data.gene_fam,
                              right_on=df_fam.gene_fam,
                              how="inner").drop(["key_0", "gene_fam_y"], axis=1)
        data_new.rename(columns={'gene_fam_x': 'gene_fam'}, inplace=True)

        data_new = data_new.sort_values(by='log10_pN/pS_mean', ascending=False)
        # to do: call figure here
        measure = "log10_pN/pS"
        self.make_box_swarm_plot(data_new, measure)

        df_top_fam = df_fam.tail(5)

        data = data.merge(df_top_fam,
                          left_on=data.gene_fam,
                          right_on=df_top_fam.gene_fam,
                          how="inner").drop(["key_0", "gene_fam_y"], axis=1)
        data.rename(columns={'gene_fam_x': 'gene_fam'}, inplace=True)

        self.plot_family_and_ref(data=data, kind="swarm", ds_order=df_top_fam.gene_fam)
        self.plot_family_and_ref(data=data, kind="box", ds_order=df_top_fam.gene_fam)

    # to do: better name
    def plot_family_and_ref(self, data, kind, ds_order):

        sns.catplot(x="gene_fam", y="log10_pN/pS", kind=kind, data=data, hue="ref",
                    palette=sns.color_palette(n_colors=10),
                    order=ds_order, height=3.5, aspect=3)

        title ="top 5 most conserved gene families {breadth}/{depth}x".format(
            breadth=self.threshold_breadth, depth=self.threshold_depth)
        plt.title(title)

        figure_name = "{}{}family_plots.families_ranked.{kind}.{title}.{depth}x.svg".format(
            self.plot_dir, self.dir_sep,
            title=title.replace(" ", "_").replace("/", "_"),
            kind=kind,
            depth=self.threshold_depth, breadth=self.threshold_breadth,
        )

        plt.savefig(figure_name)
        plt.clf()

    def quick_and_dirty_family_analysis(self):

        # now read all the files and show how the gene families vary
        # can I do this in one plot

        # *_GenePlots/crassphage_pN_pS_values.0.2.10x.txt
        # *_GenePlots/crassphage_pN_pS_values.0.8.10x.txt
        # *_GenePlots/crassphage_pN_pS_values.0.95.10x.txt

        subfolders = [f.path for f in os.scandir(self.sample_dir) if f.is_dir() and "_GenePlots" in f.name]
        samples = []

        for subfolder in subfolders:

            suffix = "{tb}.{td}x".format(tb=self.threshold_breadth, td=self.threshold_depth)
            file_name = subfolder + self.dir_sep + "crassphage_pN_pS_values.{suffix}.txt".format(suffix=suffix)

            if os.path.isfile(file_name) and os.path.getsize(file_name) > 0:
                logging.info("processing {}".format(file_name))

                temp_df = pd.read_csv(file_name,
                                      sep='\t',
                                         )
                samples.append(temp_df)

        if len(samples) == 0:
            return

        self.all_scores_df = pd.concat(samples).reset_index()

        family_df = self.all_scores_df[self.all_scores_df.gene_fam.isnull() == False]

        family_df = family_df.sort_values(by='log10_pN/pS_mean', ascending=False)
        # family_df["long_annot"] = family_df.gene_fam + family_df.Annotation

        measure = "log10_pN/pS_mean"
        self.make_box_swarm_plot(family_df, measure)

    def make_box_swarm_plot(self, family_df, measure):

        title = "Family values"
        plt.figure(figsize=(12, 10))
        plt.title("{title} ({breadth}/{depth}x)".
                  format(title=title, depth=self.threshold_depth, breadth=self.threshold_breadth))

        sns.set(style="ticks")

        agg_field = "gene_fam"
        # to do: We get a warning on the percentile calculations (implicit in box plot) for the infinite values
        # we should probably recalculate p_N/p_S with a pseudocount
        sns.boxplot(x=measure, y=agg_field, data=family_df,
                    whis="range", palette="vlag")

        sns.swarmplot(x=measure, y=agg_field, data=family_df,
                      size=2, color=".3", linewidth=0)

        self.plot_dir = self.sample_dir + self.dir_sep + "FamilyPlots"
        os.makedirs(self.plot_dir, exist_ok=True)
        figure_name = "{}{}family_plots.{}.box_plot.{measure}.{title}.{breadth}.{depth}x.svg".format(
            self.plot_dir, self.dir_sep, agg_field,
            measure=measure.replace("/", "_"), title=title.replace(" ", "_"),
            depth=self.threshold_depth, breadth=self.threshold_breadth
        )
        # plt.show()
        plt.savefig(figure_name)

        return

    def do_family_analysis(self):

        self.read_all_files_for_family()

        # to do: continue with all files for family
        self.filter_gene_sample_on_sample_and_gene_coverage(filter_genes=False)

        self.breadth_statistics()

        # to do: now we should probably filter the genes!
        data = self.filtered_gene_sample_df
        breadth_field = "breadth_{depth}x".format(depth=self.threshold_depth)
        data = data[data[breadth_field] > self.threshold_breadth]
        self.filtered_gene_sample_df = data

        self.rank_gene_families()

        self.quick_and_dirty_family_analysis()

    #endregion


def do_analysis(args_in):

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--sample_dir", dest="sample_dir",
                        help="sample directory with samples in subfolders", metavar="[sample_dir]", required=True)
    parser.add_argument("-r", "--ref", dest="ref",
                        help="reference genome id", metavar="[ref]", required=False)

    parser.add_argument("-rd", "--ref_dir", dest="ref_dir",
                        help="directory reference genomes", metavar="[ref_dir]", required=True)

    parser.add_argument("-ns", "--min_nr_samples", dest="min_nr_samples", metavar="[min_nr_samples]", type=int,
                        help="minimal nr of samples for genes")

    parser.add_argument("-td", "--threshold_depth", dest="threshold_depth", metavar="[threshold_depth]", type=int,
                        help="threshold for the depth, 1,10 or 50, will be combined with threshold for breadth")

    parser.add_argument("-tb", "--threshold_breadth", dest="threshold_breadth", metavar="[threshold_breadth]",
                        type=float,
                        help="threshold for breadth, between 0 and 1, will be combined with threshold for depth")

    parser.add_argument("-f", "--family", action="store_true", default=False,
                        help="run analysis on gene families based on the pN/pS values for individual genes",
                        required=False)

    args = parser.parse_args(args_in)

    print("Start running MakeGenePlots")
    print("sample_dir: " + args.sample_dir)
    if args.ref:
        print("reference gemome id (ref): " + args.ref)

    if not args.min_nr_samples:
        args.min_nr_samples = 3
    print("minimal nr of samples = {}".format(args.min_nr_samples))
    if not args.threshold_depth:
        args.threshold_depth = 10
    if not args.threshold_breadth:
        args.threshold_breadth = 0.95
    print("threshold depth   : {depth}".format(depth=args.threshold_depth))
    print("threshold breadth : {breadth}".format(breadth=args.threshold_breadth))

    make = MakeGenePlots(args.sample_dir, args.ref_dir, args.ref, args.threshold_depth, args.threshold_breadth)

    # also add a level of analysis that goes beyond one ref?
    # it might integrate the _gene_measures.txt created earlier
    if args.family:
        make.do_family_analysis()
    else:
        make.do_analysis(args.min_nr_samples, args.ref)


if __name__ == "__main__":
    do_analysis(sys.argv[1:])

# to do; for testing, do not use in production
# sample_dir=r"D:\17 Dutihl Lab\_tools\_pipeline\ERP005989"
# ref="crassphage_refseq"
# # ref="sib1_ms_5"
# rd = r"D:\17 Dutihl Lab\source\phages_scripts\mgx\ref_seqs"
# # do_analysis(["-d", sample_dir, "-rd", rd, "-r", ref, "-ns", "2", "-td", "10", "-tb", "0.95"])
#
# # gene family analyis (with -f option)
# do_analysis(["-d", sample_dir, "-rd", rd, "-f", "-td", "10", "-tb", "0.80"])


