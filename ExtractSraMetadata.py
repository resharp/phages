import os
import logging
import pandas as pd
from xml.dom import minidom

taxon_id = "1211417" # crAsshage

if os.name == "nt":
    sample_dir = r"D:\17 Dutihl Lab\_tools\diversitools"
else:
    sample_dir = "/hosts/linuxhome/mgx/DB/MGXDB/MGXDB"

class ExtractSraMetadata:

    sample_table_name = ""
    dir_sep = ""
    taxon_id = ""

    def __init__(self, sample_dir, taxon_id):

        if os.name == 'nt':
            self.dir_sep = "\\"
        else:
            self.dir_sep = "/"

        logging.basicConfig(filename=sample_dir + self.dir_sep + ".." + self.dir_sep + "ExtractSraMetadata.log", filemode='w', format='%(asctime)s - %(message)s',
                            level=logging.DEBUG)

        self.sample_dir = sample_dir
        self.taxon_id = taxon_id
        self.sample_table_name = sample_dir + self.dir_sep + ".." + self.dir_sep + "taxon_" + taxon_id + "_counts.txt"

    def read_metadata(self, taxon_id):

        subfolders = [f.path for f in os.scandir(self.sample_dir) if f.is_dir()]

        data = []

        for subfolder in subfolders:
            sample = os.path.basename(subfolder)
            metafile_name = subfolder + self.dir_sep + sample + "_metadata"

            if os.path.isfile(metafile_name):
                logging.debug("processing {}".format(metafile_name))

                meta = minidom.parse(metafile_name)
                items = meta.getElementsByTagName('taxon')

                taxon = [
                        [i.attributes["total_count"].value
                    ,   i.attributes["name"].value ]
                    for i in items
                    if i.attributes["tax_id"].value == taxon_id ]

                if len(taxon) == 1:
                    primary_id = meta.getElementsByTagName('PRIMARY_ID')[0].firstChild.nodeValue
                    study_title = meta.getElementsByTagName('STUDY_TITLE')[0].firstChild.nodeValue

                    data.append([sample, primary_id, taxon[0][0], taxon[0][1], study_title])

        columns = ['sample', 'primary_id', 'total_count', 'taxon_name', 'study_title']

        df = pd.DataFrame(columns=columns, data=data)

        #somehow this sorting does work but is not used in to_
        # df.sort_values(by='total_count', ascending=False, inplace=True)

        df.to_csv(path_or_buf=self.sample_table_name, sep='\t', index=False)


extract = ExtractSraMetadata(sample_dir, taxon_id)

extract.read_metadata(taxon_id)
