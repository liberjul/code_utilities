def reverse_complement(seq, rna=False):
    if rna:
        comp_dict = {"A" : "U", "U" : "A", "C" : "G", "G" : "C",
                    "N" : "N", "R" : "Y", "Y" : "R", "S" : "W",
                    "W" : "S", "K" : "M", "M" : "K", "B" : "V",
                    "V" : "B", "D" : "H", "H" : "D"}
    comp_dict = {"A" : "T", "T" : "A", "C" : "G", "G" : "C",
                "N" : "N", "R" : "Y", "Y" : "R", "S" : "W",
                "W" : "S", "K" : "M", "M" : "K", "B" : "V",
                "V" : "B", "D" : "H", "H" : "D"}
    rev_seq = seq[::-1]
    new_seq = ""
    for i in rev_seq:
        new_seq = F"{new_seq}{comp_dict[i]}"
    return new_seq

def seq_kmers(seq, k, w_pos = False):
    k_list = []
    k_dict = {}
    for i in range(0, len(seq)-k+1):
        kmer = seq[i:i+k]
        if w_pos:
            k_list.append([kmer, i, i+k])
        else:
            k_list.append(kmer)
        if kmer in k_dict:
            print(F"Non-unique kmer of length {k} found, assembly may be incorrect.")
        else:
            k_dict[kmer] = None
    return k_list

def align_w_match(seq1, seq2, kmer_list):
    for i in kmer_list:
        if i in seq1:
            kmers_w_pos1 = seq_kmers(seq1, len(i), w_pos = True)
            j = 0
            while i != kmers_w_pos1[j][0]:
                j += 1
            seq1_end = kmers_w_pos1[j][1]

            kmers_w_pos2 = seq_kmers(seq2, len(i), w_pos = True)
            j = 0
            while i != kmers_w_pos2[j][0]:
                j += 1
            seq2_start = kmers_w_pos2[j][1]
            return [seq1_end, seq2_start]

def assem_seqs(seq1, seq2, max_k = 30, min_k = 8):
    rc_seq2 = reverse_complement(seq2)
    match_found = False
    for k in range(min_k, max_k)[::-1]:
        kmers_seq2 = seq_kmers(rc_seq2, k)
        count = 0
        for i in kmers_seq2:
            if i in seq1:
                count += 1
        if count > 0:
            match_found = True
            break
    if match_found:
        match_pos = align_w_match(seq1, rc_seq2, kmers_seq2)
        assem_seq = seq1[:match_pos[0]] + rc_seq2[match_pos[1]:]
        i = 8
        while assem_seq[match_pos[0]:match_pos[0]+i] in seq1 and i < len(assem_seq):
            i += 1
        print(F"Overlap region of {i-1} bp found.")
        return assem_seq
    else:
        print("No match found with reverse complement of sequence 2. Trying without RC...")
        match_found = False
        for k in range(min_k, max_k)[::-1]:
            kmers_seq2 = seq_kmers(seq2, k)
            count = 0
            for i in kmers_seq2:
                if i in seq1:
                    count += 1
            if count > 0:
                match_found = True
                break
        if match_found:
            match_pos = align_w_match(seq1, seq2, kmers_seq2)
            assem_seq = seq1[:match_pos[0]] + seq2[match_pos[1]:]
            i = 8
            while assem_seq[match_pos[0]:match_pos[0]+i] in seq1  and i < len(assem_seq):
                i += 1
            print(F"Overlap region of {i-1} bp found.")
            return assem_seq
        else:
            assem_seq = ""
            print("No match found. Sequences do not appear to overlap with identical regions.")
            return assem_seq
