import argparse, os
import numpy as np

# Creates a consensus sequence of FASTA alignment formatted file for primer design

def consensus_code(pos_list):
    if "N" in pos_list or ("A" in pos_list and "C" in pos_list and "T" in pos_list  and "G" in pos_list):
        return "N"
    elif "A" in pos_list and "C" in pos_list and "T" in pos_list:
        return "H"
    elif "A" in pos_list and "G" in pos_list and "T" in pos_list:
        return "D"
    elif "A" in pos_list and "C" in pos_list and "G" in pos_list:
        return "V"
    elif "G" in pos_list and "C" in pos_list and "T" in pos_list:
        return "B"
    elif "A" in pos_list and "C" in pos_list:
        return "M"
    elif "A" in pos_list and "G" in pos_list:
        return "R"
    elif "A" in pos_list and "T" in pos_list:
        return "W"
    elif "C" in pos_list and "G" in pos_list:
        return "S"
    elif "C" in pos_list and "T" in pos_list:
        return "Y"
    elif "G" in pos_list and "T" in pos_list:
        return "K"
    elif "A" in pos_list:
        return "A"
    elif "C" in pos_list:
        return "C"
    elif "G" in pos_list:
        return "G"
    elif "T" in pos_list:
        return "T"
    else:
        return "-"

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="FASTA file inputs")
parser.add_argument("-o", "--output", type=str, help="FASTA file output directory")
args = parser.parse_args()

seq_dict = {}
with open(args.input, "r") as ifile:
    line = ifile.readline()
    print(line)
    while line != "":
        header = line
        line = ifile.readline()
        seq = ""
        while line != "" and line[0] != ">":
            seq += line.strip().upper()
            line = ifile.readline()
        seq_dict[header] = np.array(list(seq))

seq_arr = np.array(list(seq_dict.values()))
cons_seq = ""
for pos in range(len(seq_arr[0])):
    cons_seq += consensus_code(seq_arr[:,pos])

print(cons_seq)
with open(args.output, "w") as ofile:
    ofile.write(F">{os.path.basename(args.output).split('.')[0]}\n{cons_seq}\n")
