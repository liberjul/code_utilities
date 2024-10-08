#!/usr/bin/env python

import argparse, subprocess
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help="FASTA formatted file containing one record matching -c/--contig.")
parser.add_argument("-c", "--contig", type=str, help="A single contig record name substring, unique in file.")
parser.add_argument("-s", "--start", type=int, default = 0, help="1-indexed start location. (default: 0)")
parser.add_argument("-p", "--stop", type=int, default = -1, help="1-indexed stop location. (default: end")
parser.add_argument("-m", "--minus_strand", action="store_true", help="If the minus strand should be returned.")
args = parser.parse_args()

header = ""
count_matches = 0
with open(args.fasta, "r") as ifile:
    line = ifile.readline()
    while line != "":
        if line[0] == ">":
            if args.contig in line.strip()[1:]:
                count_matches += 1
                header = line.strip()[1:]
                line = ifile.readline()
                seq = ""
                while line != "" and line[0] != ">":
                    seq += line.strip()
                    line = ifile.readline()
        line = ifile.readline()
        
if count_matches == 0:
    raise ValueError("Contig not found in FASTA")
elif count_matches > 1:
    raise ValueError("Contig found multiple times FASTA")

if args.stop < 0:
    if args.minus_strand:
        rc_seq = str(Seq(seq[args.start-1]).reverse_complement())
        print(F">{header}_{args.start}..{len(seq)}_-\n{rc_seq}")
    else:
        print(F">{header}_{args.start}..{len(seq)}_+\n{seq[args.start-1:args.stop]}")
else:
    if args.minus_strand:
        rc_seq = str(Seq(seq[args.start-1:args.stop]).reverse_complement())
        print(F">{header}_{args.start}..{args.stop}_-\n{rc_seq}")
    else:
        print(F">{header}_{args.start}..{args.stop}_+\n{seq[args.start-1:args.stop]}")
