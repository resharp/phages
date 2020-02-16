# download metadata for samples and runs
# based on some properties
# - human gut fecal samples
# - only metagenomics data sets (and being able to filter on this in runs)
# - also add information from analysis:
#   number of reads
# e.g.
# "analysis_summary": {
# 	"Submitted nucleotide sequences": "80832",
# output: two text files
# input: please set sample_dir and lineage in script parameters
import jsonapi_client as json
import logging
import os
import pandas as pd
from urllib.parse import urlencode

API_BASE = 'https://www.ebi.ac.uk/metagenomics/api/latest/'

samples_downloaded = False

sample_dir = r"D:\17 Dutihl Lab\_tools\MGnify"

# to do change lineage
lineage = "root:Host-associated:Human:Digestive%20system:Large%20Intestine:Fecal"

if os.name == 'nt':
    dir_sep = "\\"
else:
    dir_sep = "/"

# with json.Session(API_BASE) as s:
#     study = s.get('studies', 'ERP009004').resource
#     print('Study name:', study.study_name)
#     print('Study abstract:', study.study_abstract)
#     for biome in study.biomes:
#         print('Biome:', biome.biome_name, biome.lineage)


def download_samples(lineage):

    df_samples = pd.DataFrame(columns=('sample_name', 'study', 'study_abstract'))
    df_samples.index.name = 'accession'

    df_runs = pd.DataFrame(columns=('sample_accession', 'nr_sequences','experiment_type',
                                    'instrument_platform', 'instrument model', 'analysis_pipeline'))
    df_runs.index.name = 'accession'

    with json.Session(API_BASE) as s:
        params = {
            'experiment_type': 'metagenomic'
        }
        f = json.Filter(urlencode(params))

        url = "biomes/{lineage}/samples".format(lineage=lineage)

        i_sample = 0
        for sample in s.iterate(url, f):

            # initiate sample
            logging.info("processing sample {sample}".format(sample=sample))
            i_sample = i_sample + 1

            if len(sample.studies) > 0:
                study = sample.studies[0]
                study_accession = study.accession
                study_abstract = study.study_abstract.replace("\r\n", " ").replace("\n", " ")

            df_samples.loc[sample.accession] = [
                sample.sample_name, study_accession, study_abstract
            ]
            for run in sample.runs:

                nr_sequences = 0
                for analysis in run.analyses:
                    summary = analysis.analysis_summary
                    for entry in summary:
                        if entry["key"] == "Submitted nucleotide sequences":
                            nr_sequences = entry["value"]

                df_runs.loc[run.accession] = [
                    sample.accession, nr_sequences, run.experiment_type, run.instrument_platform, run.instrument_model,
                    ", ".join([p.release_version for p in run.pipelines])
                ]

                # save some results in between?
                if i_sample % 100 == 0:
                    df_samples.to_csv(sample_dir + dir_sep + "fecal_samples.txt", sep="\t")
                    df_runs.to_csv(sample_dir + dir_sep + "fecal_runs.txt", sep="\t")

    logging.info("starting to save final data sets for samples and runs")
    df_samples.to_csv(sample_dir + dir_sep + "fecal_samples.txt", sep="\t")
    df_runs.to_csv(sample_dir + dir_sep + "fecal_runs.txt", sep="\t")

    logging.info("end of download_samples()")


logging.basicConfig(filename=sample_dir + dir_sep + 'DownloadMgnifySamples.log', filemode='w',
                    format='%(asctime)s - %(message)s', level=logging.INFO)

if not samples_downloaded:

    print("starting to download some metadata for samples/runs for lineage ${lineage}".format(lineage=lineage))
    download_samples(lineage)
else:
    print("we are NOT going to download samples")
    # just read DataFrame and get the other stuff?
