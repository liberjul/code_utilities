import argparse, glob, sys, os
import pandas as pd
from abi_handling import gen_cut_fastas_phred

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="Unzipped directory containing .ab1 and .csv files")
args = parser.parse_args()

parent = args.dir
if parent[-1] != "/":
    parent += "/"
seq_files = glob.glob(parent + "*.ab1") # Find all the chromatogram files.

# The index at which the plate position is found. For a filename like CAG384_PREMIX_CAG384_4.ab1,
# the barcode (CAG384) is at index 0 and 2 after splitting at underscores.
index_of_pos = 0 # Check that this is correct here and in the cell output

seq_dict = {}
for i in seq_files:
    seq_dict[os.path.basename(i).split("_")[index_of_pos]] = i
print(seq_dict)

meta = glob.glob(parent + "*.csv")

if len(meta) > 1:
    print("There is more than 1 metadata (*.csv) file, please remove additional.")
    sys.exit(1)
elif len(meta) == 0:
    print("There is no metadata (*.csv) file, please add it with columns 'Barcode' and 'Name'.")
    sys.exit(1)

data = pd.read_csv(meta[0])
if "Barcode" not in data.columns or "Name" not in data.columns:
    print(F"The CSV file {meta[0]} does not contain the columns 'Barcode' and 'Name'.")
    sys.exit(1)

chromat_files = []
for i in range(len(data)):
    pos = data.Barcode[i]
    os.rename(seq_dict[pos], F"{parent}{data.Name[i]}_{pos}.ab1")
    chromat_files.append(F"{parent}{data.Name[i]}_{pos}.ab1")
for i in chromat_files:
    gen_cut_fastas_phred(i, parent, os.path.basename(i).split(".")[0])
