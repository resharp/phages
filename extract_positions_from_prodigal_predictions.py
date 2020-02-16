import os
import pandas as pd

if os.name == "nt":
    mydir = r"D:\17 Dutihl Lab\crassphage\guerin"
else:
    mydir = "/hosts/linuxhome/chaperone/tmp/richard/guerin_data"

if os.name == 'nt':
    dir_sep = "\\"
else:
    dir_sep = "/"

print("extract positions from prodigal predictions")
print("processing directory: " + mydir)

#loop through all files
subfolders = [f.path for f in os.scandir(mydir) if f.is_dir()]

print("number of folders: " + str(len(subfolders)))

for subfolder in subfolders:
    # print(subfolder)
    ref = os.path.basename(subfolder)
    # print(sample)
    output_file = subfolder + dir_sep + ref + "_codingregions.txt"
    # print(output_file)

    # now get the .genes files
    files = [f.path for f in os.scandir(subfolder) if f.is_file() and "genes" in f.name]
    gene_file = ""
    if len(files) > 0:
        gene_file = files[0]
    else:
        raise NameError("No .genes file in directory " + subfolder)

    # print(gene_file)
    print("--")
    print("Starting to read {gene_file}".format(gene_file=os.path.basename(gene_file)))

    # re-start df
    df = pd.DataFrame(columns=("Beg", "End", "Reference"))
    df.index.name = "Protein"

    # extract positions and write to file

    with open(gene_file) as f:
        lines = f.readlines()
        print("read {nr_lines} lines".format(nr_lines=len(lines)))
        cds_lines = [l for l in lines if "CDS" in l]
        print("read {nr_cds} cds".format(nr_cds=len(cds_lines)))

        nr_pos = 0

        for cds in cds_lines:
            nr_pos = nr_pos + 1
            gene = ref + "_" + str(nr_pos)

            cds = cds.replace("\n", "")
            # print(cds)

            words = cds.split()
            pos_word = words[1]
            # print(pos_word)

            if pos_word.__contains__("complement"):
                pos_word = pos_word.replace("complement(", "").replace(")", "")
                positions = pos_word.split("..")
                # we reverse the positions!
                start = positions[1]
                end = positions[0]
            else:
                positions = pos_word.split("..")
                start = positions[0]
                end = positions[1]

            if ">" not in start and "<" not in start and "<" not in end and ">" not in end:
                df.loc[gene] = [start, end, ref]

    print("write {nr_cds} gene positions".format(nr_cds=len(df)))

    df.to_csv(path_or_buf=output_file, sep="\t")
