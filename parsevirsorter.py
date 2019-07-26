#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 15:35:25 2018

@author: hielke
"""

import argparse
import fileinput
import os
import re
import sys
import traceback
from collections import defaultdict
from functools import partial
from uuid import uuid4

predict_ext = "_VIRSorter_mga_final.predict"
csv_ext = "_VIRSorter_global-phage-signal.csv"

ap = argparse.ArgumentParser()
ap.add_argument(
    "inputfiles",
    nargs="*",
    help="Provide csv files here. "
    + "Predict files should be in the same directory. "
    + "CSV File is expect as {genome_id}"
    + predict_ext
    + ". "
    + "Predict File is expected as {genome_id}"
    + csv_ext
    + ".",
)
ap.add_argument("-t", "--field-seperator", default="\t")
ap.add_argument(
    "-o", "--out", nargs="?", default=sys.stdout, type=argparse.FileType("a+")
)
ap.add_argument(
    "-l", "--log", nargs="?", default=sys.stderr, type=argparse.FileType("a+")
)
args = ap.parse_args()

out = partial(print, file=args.out, sep=args.field_seperator)
log = partial(print, file=args.log, sep="\n", end="\n---- ----\n")

ID = uuid4()
log("[INFO] Program %s started" % ID)


def parse_predict_file(predict_file):
    # keys are contigs, list of all genes with a tuple(min, max)
    contig_genes = defaultdict(list)
    contig = None
    prog = re.compile(r"(\d+)")
    do_first_check = None
    with open(predict_file) as f:
        for line in f:
            if line.startswith(">"):
                found_contig = False
                while not found_contig:
                    while not line.startswith(">"):
                        line = next(f)
                    try:
                        contig = line.split("___", maxsplit=1)[0]
                        contig = contig.split("_", maxsplit=1)[1]
                        found_contig = True
                    except:
                        log(
                            "[WARING] Could not find contig in this line:",
                            line,
                            traceback.format_exc(),
                        )
                        line = next(f)

                contig_genes_list = contig_genes[contig]
                contig_genes_list.append((1, 1))
                do_first_check = True
                continue
            if contig is None:
                log("[WARNING] File did not start with >:", predict_file)
            args = line.split()
            if do_first_check is True:
                do_first_check = False
                starts_with = int(args[0].split("_")[1])
                for _ in range(starts_with - 1):
                    contig_genes_list.append((1, 1))
            if do_first_check is None:
                log("[WARNING] Something wrong with first check:", predict_file)
                continue
            if len(args) < 3:
                log("[WARNING] Something wrong with args:", predict_file)
                continue
            args = [args[1], args[2], args[-2], args[-3]]
            args = list(
                map(lambda d: int(d.group(1)), filter(bool, map(prog.match, args)))
            )
            contig_genes[contig].append((min(args), max(args)))
    return contig_genes


prog = re.compile(r"gene_(\d+)-gene_(\d+)")
for line in fileinput.input(args.inputfiles):
    if fileinput.isfirstline():
        path_name = fileinput.filename()
        genome = os.path.basename(path_name).rstrip(csv_ext)
        directory = os.path.dirname(path_name)
        predict_file = os.path.join(directory, genome + predict_ext)

        try:
            contig_genes = parse_predict_file(predict_file)
        except:
            log(
                "[FATAL] Error while parsing predict file",
                "This genome will be left unprocessed",
                predict_file,
                traceback.format_exc(),
            )
            fileinput.nextfile()
            continue

        cat = None

    if line.startswith("#"):
        c = line[3]
        if c != "C" and c in map(str, range(7)):
            cat = c
        continue
    s = sys.maxsize
    e = 0
    contig = line.split("___", maxsplit=1)[0]
    contig = contig.split("_", maxsplit=1)[1]
    for group in prog.finditer(line):
        s_pot = int(group.group(1))
        e_pot = int(group.group(2))
        if s_pot > e_pot:
            s_pot, e_pot = e_pot, s_pot
        s = min(s, s_pot)
        e = max(e, e_pot)
    if s == sys.maxsize or e == 0:
        if cat in map(str, range(4)):
            # Just the whole contig.
            out(genome, contig, cat, *4 * ("NA",))
            continue
        log("[WARNING] Cannot parse", line.strip(), genome)
        continue
    contig_genes_list = contig_genes.get(contig)
    if not contig_genes_list:
        log("[WARNING] Cannot find contig:", contig, genome)
        continue
    if s == 0:
        s = 1
    coords_contig_genes = contig_genes[contig]

    length = len(coords_contig_genes) - 1
    if e >= length:
        e = length
    else:
        e += 1

    c_s = coords_contig_genes[s][0]  # 0 is minimum
    c_e = coords_contig_genes[e][1]  # 1 is maximum

    if not cat:
        log("[WARNING] Cannot find cat:", line.trim(), genome)
        cat = "NULL"
    #        print(f'{c_s},{c_e},{s},{e}')
    out(genome, contig, cat, c_s, c_e, s, e)
log("[INFO] Program %s finished" % ID)