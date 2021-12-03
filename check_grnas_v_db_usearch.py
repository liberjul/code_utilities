#!/usr/bin/env python

# inputs:
# - sgRNA sequences, in multi-record FASTA
# - potential off-target subjects, in FASTA
# - match parameters
# - PAM sequence (NGG, for example)

# outputs
# - coordinates of mismatch hits, in gff(?) format

import glob, sys, os, argparse, subprocess, time
from auto_sanger_seq_assem import reverse_complement

def mm_and_indel(x, y):
    y = y.replace("U","T")
    total_mm = sum(c1!=c2 for c1,c2 in zip(x,y))
    indels = x.count("-") + y.count("-")
    return [total_mm-indels, indels]

def detect_strand(template, grna):
    pos = template.find(grna)
    if pos != "-1":
        return "+"
    else:
        pos = template.find(reverse_complement(grna))
        if pos != "-1":
            return "-"
        else:
            raise ValueError(F"gRNA with sequence {grna} not found in template.")
# Provide path to files and path to new, modified files
start = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--database", type=str, help="FASTA database to check for hits")
parser.add_argument("-t", "--template", type=str, help="FASTA of template to find direction of gRNAs")
parser.add_argument("-g", "--grnas", type=str, help="FASTA of gRNAs")
parser.add_argument("-u", "--usearch_path", type=str, help="Path of USEARCH executable")
parser.add_argument("-p", "--pam", default="NGG", type=str, help="PAM sequence, such as NGG")
parser.add_argument("-s", "--seed_len", default=10, type=int, help="Length of the seed sequence")
parser.add_argument("--max_seed_mismatch", default=0, type=int, help="Number of allowed seed sequence mismatches")
parser.add_argument("--max_seed_indel", default=0, type=int, help="Number of allowed seed sequence insertions/deletions")
parser.add_argument("--max_distal_mismatch", default=0, type=int, help="Number of allowed distal sequence mismatches")
parser.add_argument("--max_distal_indel", default=0, type=int, help="Number of allowed distal sequence insertions/deletions")
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

with open(args.template, "r") as ifile:
    line = ifile.readline()
    header = line.strip()[1:]
    line = ifile.readline()
    template_seq = ""
    while line != "" and line[0] != ">":
        template_seq = F"{seq}{line.strip().upper()}"
        line = ifile.readline()

PAM_str_f = args.pam
PAM_str_r = reverse_complement(args.pam)

with open(F"{args.grnas.split('.fasta')[0]}_with_PAM.fasta", "w") as ofile:
    for grna in grna_dict:
        dir = detect_strand(template_seq, grna_dict[grna])
        if dir == "+":
            seq = grna_dict[grna] + PAM_str_f
        elif dir == "-":
            seq = PAM_str_r + reverse_complement(grna_dict[grna])
        ofile.write(F">{grna}\n{seq}\n")

id_thresh = 1 - ((args.max_seed_mismatch + args.max_distal_mismatch)/max_grna_len)

print(F"Searching for matches to {len(grna_dict.keys())} gRNAs with identity >= {id_thresh}...")

usearch_comd_list = [args.usearch_path, '-usearch_global', F"{args.grnas.split('.fasta')[0]}_with_PAM.fasta", '-db', args.database, '-id', str(id_thresh), '-strand','both', '-maxaccepts','100', '-maxrejects', '0', '-userout', F"{args.output}.usearch.txt", '-userfields', 'query+target+id+mism+qstrand+tstrand+qrowdots+qrow+trow+tlor+thir']

subprocess.run(" ".join(usearch_comd_list), shell=True)

grna_hit_dict = {}
with open(F"{args.output}.usearch.txt", "r") as ifile:
    line = ifile.readline()
    while line != "":
        spl = line.strip().split("\t")
        grna_name = spl[0]
        # target_name, id, mm_count, qstrand, tstrand, query, target, tstart, tend = spl[1:]
        if grna_name not in grna_hit_dict:
            grna_hit_dict[grna_name] = [spl[1:]]
        else:
            grna_hit_dict[grna_name].append(spl[1:])
        line = ifile.readline()

# print(grna_hit_dict)

count_above_thresh = 0
oline_dict = {}
buffer ="gRNA_name,target_name,gRNA_sequence,target_sequence,mismatch_string,seed_mismatch_score,seed_indel_score,distal_mismatch_score,distal_indel_score,strand\n"#,start,stop\n"
for grna in grna_hit_dict:
    grna_hit_list = grna_hit_dict[grna]
    for hit in grna_hit_list:
        target_name, id, mm_count, qstrand, tstrand, mismatch_string, query, target, tstart, tend = hit
        if qstrand == "+":
            query_trim = query[:-len(args.pam)]
            target_trim = target[:-len(args.pam)]
            mismatch_string_trim = mismatch_string[:-len(args.pam)]
            seed_q = query_trim[-args.seed_len:]
            seed_t = target_trim[-args.seed_len:]
            distal_q = query_trim[:-args.seed_len]
            distal_t = target_trim[:-args.seed_len]
        else:
            query_trim = query[len(args.pam):]
            target_trim = target[len(args.pam):]
            mismatch_string_trim = mismatch_string[len(args.pam):]
            seed_q = query_trim[:args.seed_len]
            seed_t = target_trim[:args.seed_len]
            distal_q = query_trim[args.seed_len:]
            distal_t = target_trim[args.seed_len:]
        seed_mm, seed_indel = mm_and_indel(seed_q, seed_t)
        distal_mm, distal_indel = mm_and_indel(distal_q, distal_t)
        if seed_mm <= args.max_seed_mismatch and distal_mm <= args.max_distal_mismatch and seed_indel <= args.max_seed_indel and distal_indel <= args.max_distal_indel:
            oline = F"{grna},{target_name},{query_trim},{target_trim},{mismatch_string_trim},{seed_mm},{seed_indel},{distal_mm},{distal_indel},{qstrand}\n"
            if oline not in oline_dict:
                count_above_thresh += 1
                buffer = F"{buffer}{oline}"
                oline_dict[oline] = None

with open(F"{args.output}_hits_above_threshold.csv", "w") as ofile:
    ofile.write(buffer)

print(F"Detected {count_above_thresh} hits in {time.time()-start} seconds.")
