from Bio import SeqIO
import os

if os.name == "nt":
    mydir = r"D:\17 Dutihl Lab\crassphage"
else:
    mydir = "/hosts/linuxhome/mutant31/tmp/richard/crassphage_samples/MGXDB008660"

if os.name == 'nt':
    dir_sep = "\\"
else:
    dir_sep = "/"

gbk_file_name = mydir + dir_sep + "crassphage.gbk"

if os.path.isfile(gbk_file_name):
    parser = SeqIO.parse(gbk_file_name, "genbank")

    for rec in parser:
        for feat in rec.features:
            if feat.type == "CDS":
                locus_tag = feat.qualifiers["locus_tag"][0]
                start = feat.location.nofuzzy_start
                end = feat.location.nofuzzy_end
                product = feat.qualifiers["product"][0]

                print("{}\t{}\t{}\t{}\t{}".format(locus_tag, start, end, rec.id, product))

