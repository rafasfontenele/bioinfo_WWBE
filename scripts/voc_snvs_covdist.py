import pandas as pd
#import pickle
import sys
import os

version = sys.argv[1]

locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
mainDir = os.path.abspath(os.path.join(locDir, os.pardir))
file_path = os.path.join(mainDir, "usher_process", version, "")
data = pd.read_pickle(file_path+"lineagePaths_edited_"+version+".pkl") #, sep="\t", header=None, names=("SNV", "VOC"))
data = data.explode("SNVs")

lineages = []
ref_id_list = []
SNV = []
ref = []
alt = []
pos = []
new_data = pd.DataFrame()
metadata = pd.DataFrame()
for i, row in data.iterrows():
    if row["Lineage"].startswith("AY.") or row["Lineage"].startswith("BA.") or row["Lineage"] == "B.1.617.2" or row["Lineage"] == "B.1.1.529" \
        or row["Lineage"].startswith("BE.") or row["Lineage"].startswith("BF.") or row["Lineage"].startswith("BK.") or row["Lineage"].startswith("BJ.") \
            or row["Lineage"].startswith("BH.") or row["Lineage"].startswith("BG.") or row["Lineage"].startswith("BD.") or row["Lineage"].startswith("BC."):
        SNV.append(row["SNVs"])
        lineages.append(row["Lineage"])
        ref_id = "NC_045512.2"
        ref_id_list.append(ref_id)
        ref.append(row["SNVs"][0])
        alt.append(row["SNVs"][-1])
        pos.append(row["SNVs"][1:-1])

new_data["Ref_id"] = ref_id_list
new_data["SNV"] = SNV
new_data["Ref"] = ref
new_data["Pos"] = pos
new_data["Alt"] = alt
new_data["VOC"] = lineages

lineages_names = list(set(lineages))
metadata["names"] = lineages_names
type = []
date = []
location = []
for i, row in metadata.iterrows():
    if row['names'].startswith("AY.") or row["names"] == "B.1.617.2":
        type.append("Delta")
        location.append("Delta")
        date.append("2021-07-01")
    if row['names'].startswith("BA.") or row["names"] == "B.1.1.529" or row["names"].startswith("BE.") or row["names"].startswith("BF.") or row["names"].startswith("BK.") or row["names"].startswith("BJ.") \
            or row["names"].startswith("BH.") or row["names"].startswith("BG.") or row["names"].startswith("BD.") or row["names"].startswith("BC."):
        type.append("Omicron")
        location.append("Omicron")
        date.append("2021-12-01")

metadata["type"] = type
metadata["location"] = location
metadata["date"] = date

final_path = os.path.join(mainDir, "cov-dist_process", "")
new_data.to_csv(final_path+"pcoa_snvs_formatted_"+str(version)+".tsv", sep="\t", index=None)
metadata.to_csv(final_path+"metadata_voc_snvs_"+str(version)+".tsv", sep="\t", index=None)