#!/usr/bin/env python
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="CSV file with headers: link,gene_name. Links are to JGI's gene features, for example: \nhttps://mycocosm.jgi.doe.gov/cgi-bin/dispTranscript?db=Rhoglu1&table=protein&id=801321&useCoords=1&width=200&padding=200")
parser.add_argument("-o", "--output", type=str, help="Output CSV file for input into benchling's feature libraries")
args = parser.parse_args()

with open(args.input, "r", encoding="utf-8") as infile:
    line = infile.readline()
    line = infile.readline()
    feature_dict = {}
    browser = webdriver.Chrome()
    while line != "":
        spl_line = line.strip().split(",")
        link, alias = spl_line
        org_name = link.split("db=")[1].split("&table")[0]
        gene_id = link.split("id=")[1].split("&use")[0]
        if alias != "":
            gene_id = alias
        browser.get(link)
        red = browser.find_elements_by_xpath("/html/body/div[2]/table/tbody/tr/td[1]/form/pre/font")
        black = browser.find_elements_by_xpath("/html/body/div[2]/table/tbody/tr/td[1]/form/pre")
        whole_seq = ""
        for i in black[0].text.split("\n"):
            whole_seq += i.split(" ")[-2]
        read_dict = {}
        for i in range(0, len(red)):
            read = red[i].text
            color = red[i].get_attribute("color")
            if color == "#ff0000":
                read_dict[i] = [read, "CDS"]
            elif color == "#0000ff":
                read_dict[i] = [read, "UTR"]
            else:
                read_dict[i] = [read, "upstream"]
        # print(read_dict)
        i = 0
        while read_dict[i][1] == "upstream":
            i += 1
        utr = ""
        while read_dict[i][1] == "UTR":
            if utr + read_dict[i][0] in whole_seq:
                utr += read_dict[i][0]
            else:
                feature_dict[F"{org_name}_{gene_id}_5'UTR"] = utr
            i += 1
        if utr != "":
            feature_dict[F"{org_name}_{gene_id}_5'UTR"] = utr
        exon = ""
        exon_count = 1
        while read_dict[i][1] == "CDS":
            if exon + read_dict[i][0] in whole_seq:
                exon += read_dict[i][0]
            else:
                feature_dict[F"{org_name}_{gene_id}_exon{exon_count}"] = exon
                exon = read_dict[i][0]
                exon_count += 1
            i += 1
        feature_dict[F"{org_name}_{gene_id}_exon{exon_count}"] = exon
        utr = ""
        while read_dict[i][1] == "UTR":
            if utr + read_dict[i][0] in whole_seq:
                utr += read_dict[i][0]
            else:
                feature_dict[F"{org_name}_{gene_id}_3'UTR"] = utr
            i += 1
        if utr != "":
            feature_dict[F"{org_name}_{gene_id}_3'UTR"] = utr
        line = infile.readline()
# print(feature_dict)
with open(args.output, "w") as outfile:
    outfile.write("Name,Feature,Type\n")
    for i in feature_dict:
        if "UTR" in i:
            outfile.write(F"{i},{feature_dict[i]},UTR\n")
        else:
            outfile.write(F"{i},{feature_dict[i]},CDS\n")
