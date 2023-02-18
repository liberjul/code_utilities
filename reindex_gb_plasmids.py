#!/usr/bin/env python
import glob, sys, os, argparse
from Bio import SeqIO

'''
example usage:
python ./reindex_gb_plasmids.py \
    -i "2023_02_18_Withers_x7a_results/Withers_x7a_results/*.gbk" \
    -s ttaagattgaatcctgttgccggtcttgcgatgattatcatat # From the NOS terminator
'''

def _reindex(seq_record, loc):
    '''
    Location is zero-indexed.
    '''
    if loc > len(seq_record):
        raise ValueError("Location for reindexing is beyond the end of the sequence.")
    else:
        left = seq_record[:loc]
        right = seq_record[loc:]
        return right + left



parser = argparse.ArgumentParser()
parser.add_argument("-d", "--db", type=str, default=None, help="CSV database file for annotation records, including headers Name,Feature,Type,Color,Match type")
parser.add_argument("-i", "--input", type=str, help="Input files, wildcard can be included.")
parser.add_argument("-s", "--start_string", type=str, help="Start string, at least 10 nucleotides and expected once in a sequence.")
args = parser.parse_args()

if len(args.start_string) < 10:
    raise ValueError("The argument -s/--start_string needs to be 10+ nucleotides in length.")


files = glob.glob(args.input)
files = [x for x in files if "_reindexed.gbk" not in x]
# search_str_res = {0 : "In forward", 1: "In reverse"}
for file in files:
    res = -1
    rec = SeqIO.read(open(file, "r"), "genbank")
    id, name, description, annots, dbxrefs = rec.id, rec.name, rec.description, rec.annotations.copy(), rec.dbxrefs[:]
    if args.start_string.upper() in rec:
        res = 0
        loc = rec.seq.find(args.start_string.upper())
        new = _reindex(rec, loc)
    elif args.start_string.upper() in rec.reverse_complement():
        res = 1
        loc = rec.reverse_complement().seq.find(args.start_string.upper())
        new = _reindex(rec.reverse_complement(), loc)
    else:
        print(F"Start string not found in {file}")
    if res != -1:
        new.id, new.name, new.description, new.annotations, new.dbxrefs = id, name, description, annots, dbxrefs
        # if args.db != None:
        #
        base, ext = os.path.splitext(file)
        with open(F"{base}_reindexed.gbk", "w") as ofile:
            print(F"Reindexed {file}, length = {len(new)}")
            SeqIO.write(new, ofile, "genbank")
