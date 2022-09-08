import pandas as pd
#import sys
#import json
#import pickle
#import matplotlib.pyplot as plt
#import seaborn as sns
import os
#import numpy as np
#import matplotlib as mpl
#from ast import literal_eval
#from ete3 import Tree
#from numpy.lib.utils import info

locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
mainDir = os.path.abspath(os.path.join(locDir, os.pardir))
new_path = os.path.join(mainDir, "samples_process", "analysis", "tsvfiles", "")


# get ww samples info
def sample_info(file):# add file name as argument in future
    if os.path.isfile(new_path + str(file) + "_trimmed_union_snpEff_final.tsv"):
        tsv = pd.read_csv(new_path + str(file) + "_trimmed_union_snpEff_final.tsv", sep="\t", dtype=object)
        sample_snvs = [mut for mut in tsv['Variant']]
    else:
        print("The file " + str(file) + " is not present.")
        sample_snvs = []
        tsv = pd.DataFrame()
    return sample_snvs, tsv


# get derived snvs lineage info
def lineage_info():
    #derived_snvs = pd.read_csv("child_parent_snv_info.tsv", sep="\t")
    derived_snvs = pd.read_pickle(args.file)
    total_snvs = []
    for i, row in derived_snvs.iterrows():
        for info in row["snvs"]:
            for key, item in info.items():
                total_snvs.append(key)
    snvs = pd.read_csv(args.table, sep="\t", header=None, names=["snvs", "count"])
    snvs_to_check = snvs["snvs"].to_list()
    final_snvs = set(snvs_to_check + total_snvs)
    return final_snvs


def edit_tsv(sample, date, location, city, state, all_snvs):
    snvs, tsv = sample_info(sample)
    info = []
    for snv in snvs:
        if snv in all_snvs:
            info.append("P")
        if snv not in all_snvs:
            info.append("NP")
    tsv["Sample"] = sample
    tsv["Variant_info"] = info
    tsv["State"] = state
    tsv["City"] = city
    tsv["Date"] = date
    tsv["Location"] = location
    return tsv


def main():
    meta = pd.read_csv(args.meta, sep="\t")
    meta = meta.sort_values(by=["State", "City", "Sampler_Start_Date"])
    sample_list = [sample for sample in meta["Sample_name"]]
    all_snvs = lineage_info()
    for sample in sample_list:
        date = meta.loc[meta["Sample_name"] == sample, "Sampler_Start_Date"].iloc[0]
        location = meta.loc[meta["Sample_name"] == sample, "Location"].iloc[0]
        state = meta.loc[meta["Sample_name"] == sample, "State"].iloc[0]
        city = meta.loc[meta["Sample_name"] == sample, "City"].iloc[0]
        final_file = edit_tsv(sample, date, location, city, state, all_snvs)
        final_file.to_csv(args.outpath + "/" + str(sample) + "_trimmed_union_snpEff_final_anot.tsv", sep="\t", index=False)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Generate Heatmaps with derived SNVs info from SARS-CoV-2 and wastewater samples')
    parser.add_argument('-f', '--file', help='updated file of derived snvs', required=True, dest='file', metavar='')
    parser.add_argument('-t', '--table', help='updated snvs info nextclade', required=True, dest='table', metavar='')
    parser.add_argument('-m', '--meta', help='file containing metadata info', required=True, dest='meta', metavar='')
    parser.add_argument('-o', '--outpath', help='output folder name', required=True, dest='outpath', metavar='')
    args = parser.parse_args()

    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)
    else:
        print('Outpath already exists and will be overwritten!')
    main()
