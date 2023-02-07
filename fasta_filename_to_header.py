#!/usr/bin/env python
import glob, sys, os, argparse

# Provide path to files and path to new, modified files

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="FASTA file inputs")
parser.add_argument("-o", "--output", type=str, help="FASTA file output directory")
args = parser.parse_args()

files = glob.glob(args.input)

print(files)

for file in files:
    seq = ""
    with open(file, "r") as ifile:
        line = ifile.readline()
        print(line)
        while line != "":
            line = ifile.readline()
            seq += line.strip()
    base = os.path.basename(file)
    print(base, seq)
    with open(args.output + base, "w") as ofile:
        ofile.write(F">{base.split('.')[0]}\n")
        ofile.write(seq)
