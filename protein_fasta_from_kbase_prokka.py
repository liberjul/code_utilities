import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-b", "--genbank", type=str, help="Genbank file with amino acid sequences from Prokka")
parser.add_argument("-f", "--gff", type=str, help="GFF file with gene coordinates from Prokka")
parser.add_argument("-o", "--output", type=str, help="Output prefix of protein FASTA")

args = parser.parse_args()

gene_coord_dict = {}
with open(args.gff, "r") as ifile:
    line = ifile.readline()
    while line != "":
        if line[0] != "#":
            spl = line.split("\t")
            if spl[2] == "gene":
                start, stop = spl[3:5]
                dir = spl[6]
                ID = spl[8].split(";")[0].strip("ID=")
                if dir == "-":
                    gb_name = F"complement({start}..{stop})" # The gene name found in
                else:
                    gb_name = F"{start}..{stop}"
                gene_coord_dict[gb_name] = [start, stop, dir, ID]
        line = ifile.readline()
buffer = ""
with open(args.genbank, "r") as ifile:
        line = ifile.readline()
        while line != "" and "gene" not in line:
            line = ifile.readline()
        while line != "":
            if "gene " in line:
                line_con = ' '.join(line.split())
                print(line_con)
                gb_name = line_con.split()[1].strip()
                line = ifile.readline()
                while line != "" and ".." not in line:
                    if "product" in line:
                        prod = line.strip('"\n').split('="')[1]
                    elif "translation" in line:
                        trans = line.strip('\n').strip('"').split('="')[1]
                        line = ifile.readline()
                        while line != "" and "/" not in line and ".." not in line:
                            trans += ''.join(line.strip('\n').strip('"').split())
                            line = ifile.readline()
                    line = ifile.readline()
                buffer = F"{buffer}>{gene_coord_dict[gb_name][3]} {prod}\n{trans}\n"
            line = ifile.readline()
with open(args.output, "w") as ofile:
    ofile.write(buffer)
