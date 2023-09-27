#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--transcripts", type=str, help="path to transcripts as FASTA")
parser.add_argument("-o", "--out", type=str, help="path to output file")
args = parser.parse_args()


def seq_to_codons(seq, codon_freq):
    if len(seq) % 3 == 0:
        # print("Full length")
        for i in range(0, len(seq), 3):
            codon_freq[seq[i:i+3]] += 1
    elif seq[:3] == "ATG":
        # print("Truncated")
        for i in range(0, len(seq), 3)[:-2]:
            codon_freq[seq[i:i+3]] += 1
    else:
        pass

codon_freq = {}
for i in "AGCT":
    for j in "AGCT":
        for k in "AGCT":
            codon_freq[F"{i}{j}{k}"] = 0

with open(args.transcripts, "r") as ifile:
    line = ifile.readline()
    while line != "":
        # header = line[1:].split(" ")[0]
        line = ifile.readline()
        seq = ""
        while line != "" and line[0] != ">":
            seq += line.strip().upper()
            line = ifile.readline()
        seq_to_codons(seq, codon_freq)

total = sum([codon_freq[x] for x in codon_freq])
buffer = "Codon,Count,Frequency\n"
for i in "AGCT":
    for j in "AGCT":
        for k in "AGCT":
            codon = F"{i}{j}{k}"
            buffer += F"{codon},{codon_freq[codon]},{codon_freq[codon]/total}\n"

with open(args.out, "w") as ofile:
    ofile.write(buffer)
