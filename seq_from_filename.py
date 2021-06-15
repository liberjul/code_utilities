def seq_from_filename(file, header = False):
    seq = ""
    with open(file, "r") as ifile:
        line = ifile.readline()
        head = line
        while line != "":
            line = ifile.readline()
            seq += line.strip()
    if header:
        return [head, seq]
    else:
        return seq
