#!/usr/bin/env python

# inputs:
# - sgRNA sequences, in multi-record FASTA
# - potential off-target subjects, in FASTA
# - match parameters
# - PAM sequence (NGG, for example)

# outputs
# - coordinates of mismatch hits, in gff(?) format

# algorithm
# 1A. Create databases of all possible sgRNA target sites,
#     on both sense and antisense
#     1Ai. Create a second database that is reverse complemented
#     1Aii. Filter database to bases that are ~25 bases upstream of a PAM
# 2A. Of the possible target sites, check against provided sgRNAs
#

import glob, sys, os, argparse, subprocess, time
from auto_sanger_seq_assem import reverse_complement
import multiprocessing as mp

def look_for_hits(grna_h, grna_dict, target_dict, dir, max_grna_len, seed_len, max_seed_mismatch, max_distal_mismatch):
    below_threshold_hits = 0
    grna_hit_dict = {}
    # for grna_h in grna_dict.keys():
    seed_g = grna_dict[grna_h][-seed_len:]
    distal_g = grna_dict[grna_h][:-seed_len]
    distal_len_delta = max_grna_len - len(grna_dict[grna_h])
    grna_hit_dict[grna_h] = {}
    for target in target_dict:
        if dir == "-":
            target_r = target
            target = reverse_complement(target)
        seed_t = target[-args.seed_len:]
        distal_t = target[distal_len_delta:-args.seed_len]
        seed_mm_score = sum(c1!=c2 for c1,c2 in zip(seed_g,seed_t))
        distal_mm_score = sum(c1!=c2 for c1,c2 in zip(distal_g,distal_t))
        if seed_mm_score <= max_seed_mismatch and distal_mm_score <= max_distal_mismatch:
            if dir == "-":
                grna_hit_dict[grna_h][target_r] = [seed_mm_score, distal_mm_score, dir]
            else:
                grna_hit_dict[grna_h][target] = [seed_mm_score, distal_mm_score, dir]
            below_threshold_hits += 1
    return (below_threshold_hits, grna_hit_dict)

def find_and_report_hit_loc(grna, grna_hit_dict, verbose, database):
    buffer = ""
    for target in grna_hit_dict[grna]:
        if args.verbose:
            print(target)
        seed_mm_score, distal_mm_score, strand = grna_hit_dict[grna][target]
        # subprocess.run(F"grep -B1 '{target}' {args.database} > temp_hits.out", shell=True)
        grep_hits = subprocess.run(F"grep -B1 '{target}' {args.database}", shell=True, stdout=subprocess.PIPE).stdout.decode("utf-8").split("\n")
        # print(grep_hits)
        i = 0
        line = grep_hits[0]
        while line != "":
            header = line.strip()[1:]
            i += 1
            line = grep_hits[i]
            seq = ""
            while line != "" and line[0] != ">":
                seq = F"{seq}{line.strip().upper()}"
                i += 1
                line = grep_hits[i]
            loc = seq.find(target, 0)
            while loc != -1:
                if loc != -1 and strand == "+":
                    buffer = F"{buffer}{grna},{header},{target},{seed_mm_score},{distal_mm_score},{strand},{loc+1},{loc+1+len(target)}\n"
                elif loc != -1 and strand == "-":
                    buffer = F"{buffer}{grna},{header},{target},{seed_mm_score},{distal_mm_score},{strand},{loc+1+len(target)},{loc+1}\n"
                loc = seq.find(target, loc+1)
    return buffer

# Provide path to files and path to new, modified files
start = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--database", type=str, help="FASTA database to check for hits")
parser.add_argument("-g", "--grnas", type=str, help="FASTA of gRNAs")
parser.add_argument("-p", "--pam", default="NGG", type=str, help="PAM sequence, such as NGG")
parser.add_argument("-s", "--seed_len", default=10, type=int, help="Length of the seed sequence")
parser.add_argument("--max_seed_mismatch", default=0, type=int, help="Number of allowed seed sequence mismatches")
parser.add_argument("--max_distal_mismatch", default=0, type=int, help="Number of allowed distal sequence mismatches")
parser.add_argument("-v", "--verbose", action="store_true", help="Print target sequences and running.")
parser.add_argument("-o", "--output", type=str, help="Prefix of output files")
args = parser.parse_args()

print(F"PAM sequence 5'->3': {args.pam}")
print(F"Seed sequence length: {args.seed_len}")
print(F"Maximum mismatches in seed sequence: {args.max_seed_mismatch}")
print(F"Maximum mismatches in distal sequence: {args.max_distal_mismatch}")

grna_dict = {}
max_grna_len = 0
with open(args.grnas, "r") as ifile:
    line = ifile.readline()
    while line != "":
        header = line.strip()[1:]
        line = ifile.readline()
        seq = ""
        while line != "" and line[0] != ">":
            seq = F"{seq}{line.strip().upper()}"
            line = ifile.readline()
        if len(seq) > max_grna_len:
            max_grna_len = len(seq)
        grna_dict[header] = seq

print(F"Searching for matches to {len(grna_dict.keys())} gRNAs...")

PAM_str_f = args.pam.replace("N", ".")
PAM_str_r = reverse_complement(args.pam).replace("N", ".")
fwd_grep_str = "."*max_grna_len + PAM_str_f
rev_grep_str = PAM_str_r + "."*max_grna_len
sed_str = "."*len(args.pam)

# Make forward and reverse potential gRNA cut sites
subprocess.run(F"grep -o '{fwd_grep_str}' {args.database} | sed 's/{sed_str}$//' | sed 's/U/T/' > {args.output}_fwd_candidates.txt", shell=True)
subprocess.run(F"grep -o '{rev_grep_str}' {args.database} | sed 's/^{sed_str}//' | sed 's/U/T/' > {args.output}_rev_candidates.txt", shell=True)

fwd_dict = {}
with open(F"{args.output}_fwd_candidates.txt", "r") as ifile:
    line = ifile.readline()
    while line != "":
        fwd_dict[line.strip().upper()] = None
        line = ifile.readline()
rev_dict = {}
with open(F"{args.output}_rev_candidates.txt", "r") as ifile:
    line = ifile.readline()
    while line != "":
        rev_dict[line.strip().upper()] = None
        line = ifile.readline()

print(F"{len(fwd_dict.keys())+len(rev_dict.keys())} candidate sites adjacent to PAM sequences found...")

if mp.cpu_count() > len(grna_dict.keys()):
    pool = mp.Pool(len(grna_dict.keys()))
else:
    pool = mp.Pool(mp.cpu_count())

fwd_res = [pool.apply(look_for_hits, args=(grna_h, grna_dict, fwd_dict, "+", max_grna_len, args.seed_len, args.max_seed_mismatch, args.max_distal_mismatch)) for grna_h in grna_dict]
rev_res = [pool.apply(look_for_hits, args=(grna_h, grna_dict, rev_dict, "-", max_grna_len, args.seed_len, args.max_seed_mismatch, args.max_distal_mismatch)) for grna_h in grna_dict]

grna_hit_dict = fwd_res[0][1]
for grna_h in grna_dict:
    for i in fwd_res[1:]:
        if grna_h in grna_hit_dict and grna_h in i[1]:
            grna_hit_dict[grna_h].update(i[1][grna_h])
        elif grna_h in i[1]:
            grna_hit_dict[grna_h] = i[1][grna_h]
    for i in rev_res:
        if grna_h in grna_hit_dict and grna_h in i[1]:
            grna_hit_dict[grna_h].update(i[1][grna_h])
        elif grna_h in i[1]:
            grna_hit_dict[grna_h] = i[1][grna_h]


buffer ="gRNA_name,subject_name,target_sequence,seed_mismatch_score,distal_mismatch_score,strand,start,stop\n"

buffer_list = [pool.apply(find_and_report_hit_loc, args=(grna, grna_hit_dict, args.verbose, args.database)) for grna in grna_hit_dict]

pool.close()
buffer += "".join(buffer_list)
hit_count = buffer.count('\n')-1

print(F"Found {hit_count} targets with {args.max_seed_mismatch} or fewer seed sequence mismatched bases and {args.max_distal_mismatch} or fewer distal sequence mismatched bases...")
with open(F"{args.output}_hits_above_threshold.csv", "w") as ofile:
    ofile.write(buffer)

print(F"Finished in {time.time()-start} seconds.")
