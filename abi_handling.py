#!/usr/bin/env python
import glob, os, time
import pandas as pd
import numpy as np
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def gen_cut_fastas(in_file, path, filename, n=5):
    '''
    Inputs name of .ab1 file as string and new file name to save .fasta file to.
    Saves a .fasta file with cut sequence and changed name.
    '''
    in_data = SeqIO.read(in_file, 'abi')
    i = 0
    n = 50
    qual = False
    while qual == False:
        i += 1
        frame = list(in_data.seq[i:i+n])
        count = frame.count("N")
        if count == 0 or i > 1200:
            qual = True
            start= i
    while qual:
        i += 1
        frame = list(in_data.seq[i:i+2])
        count = frame.count("N")
        if count > 0 or i > 1200:
            qual = False
            stop = i
    cut_data = in_data[start:stop]
    cut_data.id = filename
    SeqIO.write(cut_data, path + filename + '.fasta', 'fasta')
def gen_cut_fastas_phred(in_file, path, filename, n=50, min_length=10, output=False, q_thresh=20, quals=False):
    '''
    Inputs name of .ab1 file as string and new file name to save .fasta file to.
    Saves a .fasta file with cut sequence and changed name.
    '''
    in_data = SeqIO.read(in_file, 'abi')
    seq_str, quals = in_data.seq, in_data.letter_annotations["phred_quality"]
    blanked_seq = "".join([seq_str[i] if quals[i] >= q_thresh else " " for i in range(len(seq_str))])
    seq_list = blanked_seq.split()
    longest_seq_len = 0
    cut_seq = ""
    for i in seq_list:
        if len(i) > longest_seq_len:
            longest_seq_len = len(i)
    for i in seq_list:
        if len(i) == longest_seq_len:
            cut_seq = i
            break
    cut_data = SeqRecord(Seq(cut_seq), id = filename, description="")
    start = seq_str.find(cut_data.seq)
    stop = start + longest_seq_len
    cut_quals = quals[start, stop]
    print(longest_seq_len, filename)
    if len(cut_data) >= min_length:
        SeqIO.write(cut_data, path + filename + '.fasta', 'fasta')
    if output:
        if quals:
            return filename, cut_data.seq, quals
        else:
            return filename, cut_data.seq, quals
def rename_seqs(gen_path, data_file, path, ext):
    meta_data = pd.read_csv(gen_path + data_file)
    for j in ext:
        name = gen_path + path + '*' + j
        seqs = glob.glob(name)
        print(name, seqs)
        seq_nums = list(meta_data["Sequencing_Num"])
        seq_nums.sort()
        sam_names = list(meta_data["Sample Name"])
        #sam_names.sort()
        for i in range(len(seqs)):
            if seq_nums[i] in seqs[i]:
                new_name = gen_path + meta_data["Sample Name"][i]
                if j == '.ab1':
                    gen_cut_fastas(seqs[i], gen_path, meta_data["Sample Name"][i])
                os.rename(seqs[i], new_name+j)
                print(new_name+j)
            elif sam_names[i] in seqs[i]:
                new_name = gen_path + sam_names[i]
                if j == ".ab1":
                    gen_cut_fastas(seqs[i], gen_path, sam_names[i])
                print(new_name+j)
            else:
                print(sam_names[i],  seqs[i])
                continue
def bin_nomen(result, spaces=2):
    count_spaces = 0
    i = 0
    while count_spaces < spaces:
        if result[i] == ' ':
            count_spaces += 1
        i += 1
    return result[0:i-1]
def BLAST_it(path, E_VALUE_THRESH=0.001):
    names = []
    seq_fasta = SeqIO.read(path, 'fasta')
    if len(seq_fasta.seq) > 500:
        res = NCBIWWW.qblast('blastn', 'nr', seq_fasta.seq)
        blast_records = NCBIXML.read(res)
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if 'Uncultured' in alignment.title:
                    continue
                names.append([bin_nomen(alignment.hit_def),hsp.identities/hsp.align_length])
        return names
    else:
        return ['Bad Read']
def batch_blast(file_path):
    results = {}
    seqs = glob.glob(file_path + '/*.fasta')
    for seq in seqs:
        name = os.path.basename(seq)
        t1 = time.clock()
        results[name] = BLAST_it(seq)
        t2 = time.clock()
        print(t2-t1, name, results[name])
    return results
