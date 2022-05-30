import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--reps", type=str, help="OTU representative sequence file")
parser.add_argument("-c", "--class_file", type=str, help="CONSTAX-formatted classification")
parser.add_argument("-o", "--output", type=str, help="Output file of filtered FASTAs")
parser.add_argument("-f", "--filter", type=str, help="Filter string")
args = parser.parse_args()


otu_dict = {}
with open(args.reps, "r") as ifile:
    line = ifile.readline()
    while line != "":
        header = line.strip()[1:]
        line = ifile.readline()
        seq = ""
        while line != "" and line[0] != ">":
            seq += line
            line = ifile.readline()
        otu_dict[header] = seq

# print(otu_dict)
class_df = pd.read_csv(args.class_file, sep = "\t")

fs_list = args.filter.split(";")
if len(fs_list) % 2 != 0:
    raise ValueError("Filter string must include rank and value separated by semicolons. For example\n 'Class;Microbotryomycetes'")
for i in fs_list[::2]:
    if i not in class_df.columns:
        raise ValueError(F"Filter string ranks must be columns of the classification file: {class_df.columns[1:]}")

sub_df = class_df
for i,j in zip(fs_list[::2], fs_list[1::2]):
    if j == "":
        raise ValueError("Value critereon must not be empty. The should be a taxon name or 'NaN' preceded by nothing or '!'")
    elif j == "NaN":
        sub_df = sub_df[pd.isnull(sub_df[i])]
    elif j[0] == "!":
        j = j.strip("!").strip()
        if j != "NaN":
            sub_df = sub_df[sub_df[i] != j]
        else:
            sub_df = sub_df[pd.notnull(sub_df[i])]
    else:
        sub_df = sub_df[sub_df[i] == j]

print(sub_df)

buffer = ""
for i in sub_df["OTU_ID"]:
    buffer += F">{i}\n{otu_dict[i]}"

with open(args.output, "w") as ofile:
    ofile.write(buffer)
