#!/usr/bin/env python
import numpy as np
import argparse, time

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="FASTQ file input")
args = parser.parse_args()

decode_dict = {'!' : 0,
    '"' : 1,
    '#' : 2,
    '$' : 3,
    '%' : 4,
    '&' : 5,
    "'" : 6,
    '(' : 7,
    ')' : 8,
    '*' : 9,
    '+' : 10,
    ',' : 11,
    '-' : 12,
    '.' : 13,
    '/' : 14,
    '0' : 15,
    '1' : 16,
    '2' : 17,
    '3' : 18,
    '4' : 19,
    '5' : 20,
    '6' : 21,
    '7' : 22,
    '8' : 23,
    '9' : 24,
    ':' : 25,
    ';' : 26,
    '<' : 27,
    '=' : 28,
    '>' : 29,
    '?' : 30,
    '@' : 31,
    'A' : 32,
    'B' : 33,
    'C' : 34,
    'D' : 35,
    'E' : 36,
    'F' : 37,
    'G' : 38,
    'H' : 39,
    'I' : 40}

def decode_char(char):
    return decode_dict[char]
decode_char_vec = np.vectorize(decode_char)

def decode_q(q_line):
    line_arr = np.array(list(line))
    return decode_char_vec(line_arr)

def arr_mean(old_arr, new_arr):
    return (old_arr + new_arr)/2



with open(args.input, "r") as ifile:
    start = time.time()
    char_count = 0
    run_sum = 0
    line = ifile.readline()
    while line != "":
        for i in range(3):
            line = ifile.readline()
        line = line.strip()
        run_sum += (np.sum(decode_q(line)))
        char_count += len(line)
        line = ifile.readline()
    print(args.input, " ", run_sum/char_count)
    print("Found in ", {time.time()-start}, " seconds")
