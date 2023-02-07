#!/usr/bin/env python
import glob, sys, os, argparse

# Provide path to files and path to new, modified files

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="Directory containing FASTA files")
parser.add_argument("-i", "--ignore_order", action="store_true", help="Ignore record order (improves speed for large files).")
args = parser.parse_args()

if args.dir[-1] != "/":
    args.dir += "/"

files = glob.glob(args.dir + "*.fasta")

print(files)

for file in files:
    seq_dict = {}
    if not args.ignore_order:
        header_list = []
    with open(file, "r") as ifile:
        line = ifile.readline()
        while line != "":
            header = line
            if not args.ignore_order:
                header_list.append(header)
            seq = ""
            line = ifile.readline()
            while line != "" and line[0] != ">":
                seq += line.strip()
                line = ifile.readline()
            seq_dict[header] = seq
            line = ifile.readline()
    with open(file, "w") as ofile:
        if not args.ignore_order:
            for head in header_list:
                ofile.write(head)
                ofile.write(seq_dict[head])
        else:
            for head in seq_dict:
                ofile.write(head)
                ofile.write(seq_dict[head])
