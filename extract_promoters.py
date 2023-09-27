#!/usr/bin/env python
import argparse
from Bio.Seq import Seq

'''
/mnt/c/Users/julia/Anaconda3/Lib/site-packages/code_utilities/extract_promoters.py \
  -f JL201_WT_47C1-T2-3_flye_racon_assembly.fasta \
  -g Microbotryomycetes_sp_JL201_rRNA_added_manual_v2.gff3 \
  -o Microbotryomycetes_sp_JL201_rRNA_added_manual_v2.promoters.fna \
  -l 1000
'''

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help="FASTA file of genome")
parser.add_argument("-g", "--gff", type=str, help="GFF file with gene coordinates")
parser.add_argument("-o", "--output", type=str, help="Output FASTA of promoters")
parser.add_argument("-l", "--length", type=int, default = 1000, help="Promoter length in bp")
args = parser.parse_args()

gene_coord_dict = {}
with open(args.gff, "r") as ifile:
    line = ifile.readline()
    while line != "":
        if line[0] != "#":
            spl = line.split("\t")
            if spl[2] == "gene":
                contig = spl[0]
                start, stop = spl[3:5]
                dir = spl[6]
                ID = spl[8].split(";")[0].strip("ID=")
                # if dir == "-":
                #     gb_name = F"complement({start}..{stop})" # The gene name found in
                # else:
                #     gb_name = F"{start}..{stop}"
                gene_coord_dict[ID] = [int(start), int(stop), dir, contig]
        line = ifile.readline()

seq_dict = {}
with open(args.fasta, "r") as ifile:
    line = ifile.readline()
    while line != "":
        header = line[1:].split(" ")[0]
        line = ifile.readline()
        seq = ""
        while line != "" and line[0] != ">":
            seq += line.strip().upper()
            line = ifile.readline()
        seq_dict[header] = seq
seq_lengths = {k : len(seq_dict[k]) for k in seq_dict}

seq_IDs = [x for x in gene_coord_dict]
seq_IDs.sort()

with open(args.output, "w") as ofile:
    for i in seq_IDs:
        start, stop, dir, contig = gene_coord_dict[i] # the start and stop are 1-indexed
        if dir == "+":
            prom_start = max(1, start-args.length)
            prom_end = start - 1
        elif dir == "-":
            prom_start = stop + 1
            prom_end = min(stop+args.length, seq_lengths[contig])
        seq = seq_dict[contig][prom_start - 1:prom_end]
        if dir == "-":
            seq = str(Seq(seq).reverse_complement())
        ofile.write(F">{i}_promoter_{args.length}bp\n{seq}\n")
