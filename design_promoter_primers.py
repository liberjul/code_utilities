#!/usr/bin/env python
import primer3, argparse
import pandas as pd
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help="Promoters FASTA")
parser.add_argument("-d", "--dtable", type=str, default = "", help="Text file with expression values for each gene")
parser.add_argument("-s", "--sites", type=str, default = "GGTCTC", help="List of comma-separated restriction sites to exclude from sequences")
parser.add_argument("-o", "--output", type=str, help="Text file output")
parser.add_argument("--fwd_prefix", type=str, default = "ACGTGTGGTCTCGatct", help="5'-3' forward primer bases before binding site")
parser.add_argument("--rev_prefix", type=str, default = "ACGTGTGGTCTCGattg", help="5'-3' reverse primer bases before binding site")
args = parser.parse_args()

BsaI_sites = ["ggtctc", "gagacc"]
BsmBI_sites = ["gctctc", "gagacg"]

prom_dict = {}
with open(args.fasta, "r") as in_prom:
    line = in_prom.readline()
    while line != "":
        header = line.strip()[1:]
        line = in_prom.readline()
        seq = ""
        while line != "" and line[0] != ">":
            seq = F"{seq}{line.strip()}"
            line = in_prom.readline()
        prom_dict[header] = seq.lower()

primer_dict = {}
count = 0
for i in prom_dict:
    # if count > 100:
    #     break
    count += 1
    primer_dict[i] = []
    seq = prom_dict[i]
    while "ttg" in seq:
        prim_len = 20
        tm_fwd = 0
        while prim_len < 30 and tm_fwd < 60:
            cand_primer_fwd = seq[-prim_len:-3]
            tm_fwd = primer3.calc_tm(cand_primer_fwd, dntp_conc=0.8, dna_conc=50, dv_conc=1.5)
            # print(cand_primer_fwd, tm_fwd)
            if tm_fwd < 50:
                prim_len += 1
            else:
                cand_primer_fwd = str(Seq(cand_primer_fwd).reverse_complement())
                break
        # print(cand_primer_fwd, tm_fwd)
        dg_homo_fwd = primer3.calc_homodimer(cand_primer_fwd, dntp_conc=0.8, dna_conc=50, dv_conc=1.5).dg
        if dg_homo_fwd < -10000:
            break
        loc = seq.find("ttg")
        if len(seq) - loc > 400:
            # print("Fragment length:", len(seq) - loc)
            prim_len = 18
            tm_rev = 0
            while prim_len < 27 and tm_rev < 60:
                cand_primer_rev = seq[loc:loc+prim_len]
                tm_rev = primer3.calc_tm(cand_primer_rev, dntp_conc=0.8, dna_conc=50, dv_conc=1.5)
                if tm_rev < 50:
                    prim_len += 1
                else:
                    break
            dg_homo_rev = primer3.calc_homodimer(cand_primer_rev, dntp_conc=0.8, dna_conc=50, dv_conc=1.5).dg
            dg_hetero = primer3.calc_heterodimer(cand_primer_rev, cand_primer_fwd, dntp_conc=0.8, dna_conc=50, dv_conc=1.5).dg
            if dg_homo_rev < -10000 or dg_hetero < -10000 or abs(tm_fwd - tm_rev) > 5:
                seq = seq[loc+3:]
            else:
                # print(F"{i}: {cand_primer_fwd},{cand_primer_rev},{tm_fwd},{tm_rev},{len(seq) - loc}")
                sub_seq = seq[loc:]
                n_bsai_sites = sum([sub_seq.count(x) for x in BsaI_sites])
                n_bsmbi_sites = sum([sub_seq.count(x) for x in BsmBI_sites])
                primer_dict[i].append(F"{args.fwd_prefix+cand_primer_fwd}\t{args.rev_prefix+cand_primer_rev}\t{tm_fwd}\t{tm_rev}\t{len(seq) - loc}\t{seq[loc:]}\t{n_bsai_sites}\t{n_bsmbi_sites}")
                seq = seq[loc+3:]
        else:
            seq = ""

if args.dtable != "":
    expr_dat = pd.read_table(args.dtable)

    buffer = "promoter_header\tfwd_primer\trev_primer\tfwd_tm\trev_tm\tprod_length\tseq\tN_BsaI_Sites\tFPKM\tGoNames\n"

    for i in primer_dict:
        if len(primer_dict[i]) > 0:
            header_split = i.split("__")
            protID = int(header_split[5])
            sub_dat = expr_dat[(expr_dat.protein_id == protID)]
            for j in primer_dict[i]:
                if len(sub_dat) > 0 and "structural constituent of ribosome" not in str(sub_dat.GoNames.iloc[0]):
                    # print(sub_dat, sub_dat.protein_id.iloc[0], protID)
                    buffer += F"{i}\t{j}\t{sub_dat.FPKM.iloc[0]}\t{sub_dat.GoNames.iloc[0]}\n"
else:
    buffer = "promoter_header\tfwd_primer\trev_primer\tfwd_tm\trev_tm\tprod_length\tseq\tN_BsaI_Sites\tN_BsmBI_Sites\n"

    for i in primer_dict:
        if len(primer_dict[i]) > 0:
            header_split = i.split("__")
            for j in primer_dict[i]:
                buffer += F"{i}\t{j}\n"

with open(args.output, "w") as ofile:
    ofile.write(buffer)

print(F"{count} promoters processed...")
