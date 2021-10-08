import argparse

def build_bibtex_dict(ifile, db_dict):
    line = ifile.readline()
    while line != "":
        if line[0] == "@":
            entry_dict = {}
            entry_dict["entry_name"] = line
            line = ifile.readline()
            while line != "" and line[0] != "@":
                spl = line.split(" = ")
                if len(spl) == 2:
                    field, data = spl
                    entry_dict[field.strip()] = data
                line = ifile.readline()
            if "title" in entry_dict:
                title = entry_dict["title"].strip(",\n").strip("{").strip("{").strip("}").strip("}")
                db_dict[title] = entry_dict
            line = ifile.readline()
        else:
            line = ifile.readline()
    return db_dict

def make_bibtex_from_dict(db_dict, key_list):
    buffer = ""
    end = "}\n"
    for key in key_list:
        entry_dict = db_dict[key]
        buffer = F"{buffer}{entry_dict['entry_name']}"
        for i in entry_dict:
            if i != "entry_name":
                buffer = F"{buffer}{i} = {entry_dict[i]}"
        buffer += end
    return buffer


parser = argparse.ArgumentParser()
parser.add_argument("--db1", type=str, help="Database to take titles from and add to 2")
parser.add_argument("--db2", type=str, help="Database to add titles in 1 but not already in 2")
parser.add_argument("-o", "--output", type=str, help="Output file")
args = parser.parse_args()

db1_dict = {}
with open(args.db1, "r") as ifile:
    db1_dict = build_bibtex_dict(ifile, db1_dict)
db2_dict = {}
with open(args.db2, "r") as ifile:
     db2_dict = build_bibtex_dict(ifile, db2_dict)

keys_to_add = list(set(list(db1_dict.keys())) - set(list(db2_dict.keys())))

buffer = make_bibtex_from_dict(db1_dict, keys_to_add)

with open(args.output, "w") as ofile:
    ofile.write(buffer)
