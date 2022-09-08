#!/usr/bin/env python
#from numpy.lib.utils import info
import pandas as pd
from ete3 import Tree
#import pickle
import sys
import json
import os

# date of current version
version = sys.argv[1]
# name of the file coming from gisaid with sequence insertions and deletions "Full_table_info_ins_del_final_gtr75_listsnvs.tsv"
ins_del = sys.argv[2]
# usher processes info of snvs per lineage "lineagePaths_edited_<date>.pkl"
snvs = sys.argv[3]

locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
mainDir = os.path.abspath(os.path.join(locDir, os.pardir))
alias_path = os.path.join(mainDir, "derived_snvs", "alias_key.json")

# Load the alias info
with open(alias_path) as f:
    d = json.load(f)
    m = {key: value if value else key for key, value in d.items() if not type(value) is list}
f.close()


def get_pair(i):
    l = i.split(".")
    pair = set()
    while len(l) > 1:
        suf = l.pop()
        pre = ".".join(l)
        if len(l) == 1:
            pair.add((m[pre], pre+"."+suf))
        else:
            pair.add((pre, pre+"."+suf))
    return pair


# Load parent children key from alias_key
lineage_path = os.path.join(mainDir, "derived_snvs", "lineage_list.txt")
data = [line.strip() for line in open(lineage_path, 'r')]
t = {("root", "A"), ("root", "B")}
for i in data:
    t.update(get_pair(i))


tree = Tree.from_parent_child_table(t)

# merge snvs from Usher and ins/deletion derived from GISAID
file_ins_del = pd.read_csv(ins_del, sep="\t", header=None, names=["Lineage", "root_id", "SNVs"])
file_ins_del.set_index("Lineage", inplace=True)
file_snvs = pd.read_pickle(snvs)
file_snvs.set_index("Lineage", inplace=True)
for i in file_snvs.index:
    for index in file_ins_del.index:
        if i == index:
            list_to_add = file_ins_del.loc[index, "SNVs"].split(",")
            old_list = file_snvs.loc[i, "SNVs"]
            new_list = old_list + list_to_add
            file_snvs.at[i, "SNVs"] = new_list


# function for store info in a dataframe
def store_info(file_snvs, tree):
    data = pd.DataFrame()
    lineages_list = []
    parents_list = []
    child_list = []
    snvs_list = []
    derived_snvs_list = []
    for i, row in file_snvs.iterrows():
        sample = i
        sample_snvs = row["SNVs"]
        sample_snvs = list(set(sample_snvs))
        if sample not in lineages_list:
            children, parent, sisters = get_relationship(sample, tree)
            lineages_list.append(sample)
            parents_list.append(parent)
            child_list.append(children)
            snvs_list.append(snvs_lineage_all(sample_snvs, file_snvs))  # changed to get info of snvs shared by lineage
            list_snvs = get_defining_snvs(sample, parent, file_snvs)  # add 
            derived_snvs_list.append(snvs_lineage_sisters(list_snvs, sisters, file_snvs, sample))
    data["lineage"] = lineages_list
    data['parent'] = parents_list
    data['child'] = child_list
    data["snvs"] = snvs_list
    data["derived_snvs"] = derived_snvs_list
    return data


def check_missing_parent(file_snvs, tree):
    temp_df = pd.DataFrame()
    lineages_to_check = [lin for lin in file_snvs.index]
    to_add = []
    for i, row in file_snvs.iterrows():
        lineage = i #row["Lineage"]
        child, parent, sister = get_relationship(lineage, tree)
        if parent not in lineages_to_check and parent != "root" and parent != "B" and parent != "":
            to_add.append(parent)
    to_add = list(set(to_add))
    if to_add:
        for lin in to_add:
            snvs_list = []
            child, parent, sister = get_relationship(lin, tree)
            child_check = [ch for ch in child if ch in lineages_to_check]
            child_num = len(child_check)
            for children in child_check:
                if children in file_snvs.index:
                    info = file_snvs.loc[children, "SNVs"]
                    snvs_list += info
            snvs_list_shared= [snv for snv in snvs_list if snvs_list.count(snv) == child_num]
            temp_df = temp_df.append({"Lineage": lin, "number": 0, "SNVs": snvs_list_shared}, ignore_index=True)
    temp_df = temp_df.set_index("Lineage")
    file_snvs = file_snvs.append(temp_df)
    return file_snvs


# function to get relationship info
def get_relationship(lineage, tree):
    child = []
    sisters = []
    parent = ""
    for node in tree.traverse():
        if node.name == lineage:
            c = node.get_children()
            parent = (node.up).name
            s = node.get_sisters()
            for kid in c:
                child.append(kid.name)
            for sis in s:
                sisters.append(sis.name)
    return child, parent, sisters


# function to get derived SNVs
def get_defining_snvs(lineage, parent, file_snvs):
    parent_snvs_list = []  # for when a parent is missing from the data
    for i, row in file_snvs.iterrows():
        if i == lineage:
            lineage_snvs_list = row["SNVs"]
            lineage_snvs_list = list(set(lineage_snvs_list))
        if i == parent:
            parent_snvs_list = row["SNVs"]
            parent_snvs_list = list(set(parent_snvs_list))
    final_snvs_list = [snv for snv in lineage_snvs_list if snv not in parent_snvs_list]
    return final_snvs_list


# function to get lineages associated with derived snvs
def snvs_lineage_sisters(list_snvs, sisters, file_snvs, sample):
    snv_lin_list = []
    for snvs in list_snvs:
        info = {}
        info[snvs] = []
        for i, row in file_snvs.iterrows():
            if snvs in row["SNVs"]:
                if i in sisters or i == sample:
                    info[snvs].append(i)
        snv_lin_list.append(info)
    return snv_lin_list


def snvs_lineage_all(list_snvs, file_snvs):
    snv_lin_list = []
    for snvs in list_snvs:
        info = {}
        info[snvs] = []
        for i, row in file_snvs.iterrows():
            if snvs in row["SNVs"]:
                info[snvs].append(i)
        snv_lin_list.append(info)
    return snv_lin_list


checked_file = check_missing_parent(file_snvs, tree)
result = store_info(checked_file, tree)
final_path = os.path.join(mainDir, "derived_snvs")
result.to_pickle(final_path+"/child_parent_info_"+version+"_all.pkl")  # to load directly into lineage_assingment
result.to_csv(final_path+"/child_parent_info_"+version+"_all.tsv", sep="\t", index=False)# output_file argument 