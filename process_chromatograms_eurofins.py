#!/usr/bin/env python
import argparse, glob, sys, os
import pandas as pd
from abi_handling import gen_cut_fastas_phred

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="Unzipped directory containing .ab1 and .csv files")
parser.add_argument("-n", "--nucs_in_window", type=int, default=50, help="Nucleotides in sliding window; first stretch of at least n nucleotides with quality above 20 are output.")
parser.add_argument("-m", "--min_length", type=int, default=50, help="Minimum length of sequence to be output as a FASTA.")
parser.add_argument("-i", "--index", type=int, default=0, help="Index of barcode after splitting on underscores.")
parser.add_argument("-r", "--rename", action="store_true", default=False, help="If files should be renamed.")
args = parser.parse_args()

parent = args.dir
if parent[-1] != "/":
    parent += "/"
seq_files = glob.glob(parent + "*.ab1") # Find all the chromatogram files.

# The index at which the plate position is found. For a filename like CAG384_PREMIX_CAG384_4.ab1,
# the barcode (CAG384) is at index 0 and 2 after splitting at underscores.
# index_of_pos = 0 # Check that this is correct here and in the cell output

seq_dict = {}
for i in seq_files:
    seq_dict[os.path.basename(i).split("_")[args.index].split(".ab1")[0]] = i
print(seq_dict)

meta = [x for x in glob.glob(parent + "*.csv") if x != parent + "sequences.csv"]

no_meta  = False
if len(meta) > 1:
    print("There is more than 1 metadata (*.csv) file, please remove additional.")
    sys.exit(1)
elif len(meta) == 0:
    print("There is no metadata (*.csv) file. If you wish to rename sequences, please add a CSV file with columns 'Barcode' and 'Name'.")
    cont = input("Do you wish to continue? y/n \n").lower().strip()
    if cont == "y":
        print("Continuing...")
        no_meta = True
    else:
        print("Exiting...")
        sys.exit(1)

if no_meta:
    buffer = "Name,Sequence\n"
    for i in seq_files:
        name, seq = gen_cut_fastas_phred(i, parent, os.path.basename(i).split(".")[0], output = True, n=args.nucs_in_window, min_length = args.min_length)
        if len(seq) >= args.min_length:
            buffer = F"{buffer}{name},{seq}\n"
else:
    data = pd.read_csv(meta[0])
    if "Barcode" not in data.columns or "Name" not in data.columns:
        print(F"The CSV file {meta[0]} does not contain the columns 'Barcode' and 'Name'.")
        sys.exit(1)

    chromat_files = []
    if args.rename:
        _ = [seq_dict[x] for x in data.Barcode]
    for i in range(len(data)):
        if args.rename:
            pos = data.Barcode[i]
            os.rename(seq_dict[pos], F"{parent}{data.Name[i]}_{pos}.ab1")
            chromat_files.append(F"{parent}{data.Name[i]}_{pos}.ab1")
        else:
            chromat_files = glob.glob(parent + "*.ab1")
    buffer = "Name,Sequence\n"
    for i in chromat_files:
        name, seq = gen_cut_fastas_phred(i, parent, os.path.basename(i).split(".")[0], output = True, n=args.nucs_in_window, min_length = args.min_length)
        if len(seq) >= args.min_length:
            buffer = F"{buffer}{name},{seq}\n"

with open(parent + "sequences.csv", "w") as ofile:
    ofile.write(buffer)
