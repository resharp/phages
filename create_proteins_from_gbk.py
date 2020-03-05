from Bio import SeqIO
import os

if os.name == "nt":
    seq_dir = r"D:\17 Dutihl Lab\source\phages_scripts\mgx\ref_seqs"
    dir_sep = "\\"
else:
    seq_dir = "scripts/mgx/ref_seqs"
    dir_sep = "/"

ref = "crassphage_refseq"

gbk_filename = seq_dir + dir_sep + "crassphage.gbk".format(ref=ref)
faa_filename = seq_dir + dir_sep + "{ref}.proteins.faa".format(ref=ref)

print("create_proteins_from_gbk.py: start converting {gbk} to {faa}".format(
    gbk=gbk_filename, faa=faa_filename
))

gbk_input = open(gbk_filename, "r")
faa_output = open(faa_filename, "w")

for seq_record in SeqIO.parse(gbk_input, "genbank"):

    print("Dealing with GenBank record %s" % seq_record.id)

    for seq_feature in seq_record.features:

        if seq_feature.type == "CDS":

            assert len(seq_feature.qualifiers['translation']) == 1

            faa_output.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))

faa_output.close()

gbk_input.close()