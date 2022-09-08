import pandas as pd
import sys

a = pd.read_csv(sys.argv[1], sep="\t", header=0, names=["Lineage", "Count"])
b = pd.read_csv(sys.argv[2], sep="\t", header=0, names=["SNV", "Lineage", "Count"])

info = {}
seqs = []
prevalence = []
for i, lin in a.iterrows():
    info[lin["Lineage"]] = int(lin["Count"])

for i, line in b.iterrows():
    if line["Lineage"] in info.keys():
        seqs.append(info[line["Lineage"]])
        calc = (int(line["Count"])/info[line["Lineage"]])*100
        prevalence.append(calc)
    elif line["Lineage"] not in info.keys():
        print("Something went wrong!")

b["Total_seqs"] = seqs
b["Prevalence"] = prevalence


out = b.to_csv(sys.argv[3], index=False, sep="\t")