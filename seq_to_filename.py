def seq_to_filename(header, seq, file, wrap = False, wrap_count = 60):
    with open(file, "w") as ofile:
        ofile.write(">" + header.strip(">").strip() + "\n")
        seq = seq.strip()
        if wrap:
            while len(seq) >= wrap_count:
                ofile.write(seq[0:wrap_count] + "\n")
                seq = seq[wrap_count:]
            ofile.write(seq + "\n")
        else:
            ofile.write(F"{seq}\n")
