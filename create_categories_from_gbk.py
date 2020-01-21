from Bio import SeqIO
import os
import pandas as pd

if os.name == "nt":
    mydir = r"D:\17 Dutihl Lab\crassphage"
else:
    mydir = "/hosts/linuxhome/mutant31/tmp/richard/crassphage_samples/MGXDB008660"

if os.name == 'nt':
    dir_sep = "\\"
else:
    dir_sep = "/"

gbk_file_name = mydir + dir_sep + "crassphage.gbk"

df = pd.DataFrame(columns=("Beg", "End", "Reference"))
df.index.name = "Protein"

if os.path.isfile(gbk_file_name):
    parser = SeqIO.parse(gbk_file_name, "genbank")

    for rec in parser:
        strand = 0
        for feat in rec.features:
            if feat.type == "gene":
                # this may contain "complement([start]..[end])
                debug = True
                strand = feat.strand
            if feat.type == "CDS":
                locus_tag = feat.qualifiers["locus_tag"][0]
                start = feat.location.nofuzzy_start + 1 # +1 because the numbering in biopython seems to be 0 based
                end = feat.location.nofuzzy_end
                product = feat.qualifiers["product"][0]

                if locus_tag != "KP06_gp34":
                    if strand == 1:
                        df.loc[locus_tag] = [
                            start, end, rec.id
                        ]
                    else:
                        df.loc[locus_tag] = [
                            end, start, rec.id
                        ]
                else:
                    part_1 = locus_tag + "_1"
                    part_2 = locus_tag + "_2"
                    df.loc[part_1] = [
                        start, 24472, rec.id
                    ]
                    df.loc[part_2] = [
                        25489, end, rec.id
                    ]

# change the positions of the first gene KP06_gp01 to 191 and 346
df.loc["KP06_gp01", "Beg"] = 191
df.loc["KP06_gp01", "End"] = 346

# split KP06_gp34 in two like this:
# KP06_gp34_1	24167	24472	NC_024711.1
# KP06_gp34_2	25489	25851	NC_024711.1

# to do write to file crassphage_codingregions.txt
file_name = mydir + dir_sep + "crassphage_codingregions.txt"
df.to_csv(path_or_buf=file_name, sep="\t")


