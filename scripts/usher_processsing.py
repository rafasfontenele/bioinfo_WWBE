import numpy as np
import pandas as pd
import sys
import pickle
import subprocess
import urllib.request
import os
from ete3 import PhyloTree
import logging, time

#this function updates the Usher tree file
def download_tree(locDir):
	url = "http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/"\
	          "UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
	treePath = os.path.join(locDir, "public-latest.all.masked.pb.gz")
	urllib.request.urlretrieve(url, treePath)
	return treePath

# this function generates the files to from the Usher tree
def generate_files_from_tree(TreePath):
    cmd_lineagepath = f"matUtils extract -i {TreePath} -C lineagePaths.txt"
    cmd_allpaths = f"matUtils extract -i {TreePath} -A all_Paths.txt"
    cmd_nwktree = f"matUtils extract -i {TreePath} -t Usher_fulltree.nwk"
    sys.stdout.flush()  # force python to flush
    completed_lin = subprocess.run(cmd_lineagepath, shell=True, executable="/bin/bash",
                               stdout=subprocess.DEVNULL)
    completed_all = subprocess.run(cmd_allpaths, shell=True, executable="/bin/bash",
                               stdout=subprocess.DEVNULL)
    completed_tree = subprocess.run(cmd_nwktree, shell=True, executable="/bin/bash",
                               stdout=subprocess.DEVNULL)
    return completed_lin, completed_all, completed_tree

# this function edits the all_paths.txt 
def edit_all_paths():
    df = pd.DataFrame()
    with open("all_Paths.txt", "r") as file:
        keep=[line for line in file]
        nodes_name=[info.split(":")[0] for info in keep]
        nodes_snvs=[info.split(":")[1].strip(" ").rstrip() for info in keep]
        df["Nodes"] = nodes_name
        df["SNVs"] = nodes_snvs
    df.to_csv("node_snvs_all.tsv", sep="\t", index=False)

#auxiliar function
def fix(x):
    to_add =[]
    to_remove =[]
    for mut in x:
        if "," in mut:
            to_remove.append(mut)
            for nmut in mut.split(","):
                to_add.append(nmut)
    x = x + to_add
    x = [item for item in x if item not in to_remove]
    return x

# this function formats the lineagePaths.txt file to turn the snvs info into a list
def parse_tree_paths(df, clade):
    df = df.set_index('clade')
    #take out nextrain clade name
    if clade == "clade_out":
        nxNames = df.index[df.index.str[0].str.isdigit()]
        df = df.drop(index=nxNames)
    # Make sure to check with new tree versions, lineages could get trimmed.
    df = df.drop_duplicates(keep='last')
    df['from_tree_root'] = df['from_tree_root'].fillna('')
    df['from_tree_root'] = df['from_tree_root']\
        .apply(lambda x: x.replace(' ', '').strip('>').split('>'))
    df['from_tree_root'] = df['from_tree_root']\
        .apply(fix)
    df.reset_index(inplace=True)
    df = df.rename(columns={"clade":"Lineage", "from_tree_root":"SNVs"})
    return df

#This function (written by Yiyan) prepares the Usher tree to keep lineage node only
def create_lineage_tree_for_itol():
    cmd_info= f"cut -f2 lineagePaths.txt | tail -n +2 | sort | uniq > lineage_nodes"
    sys.stdout.flush()
    completed = subprocess.run(cmd_info, shell=True, executable="/bin/bash",
                               stdout=subprocess.DEVNULL)
    tree = PhyloTree("Usher_fulltree.nwk", format=3)
    with open("lineage_nodes") as f:
        node_list = f.read().splitlines()
    for node in tree.traverse("postorder"):
        if node.is_leaf() and node.name not in node_list:
            node.detach()
    for node in tree.traverse("postorder"):
        if node.up:
            parent = node.up
            if len(parent.get_children())==1 and (not node.name in node_list) and (not node.is_leaf()):
                for child in node.get_children():
                    parent.add_child(child)
                parent.remove_child(node)
    tree.write(format=3, outfile="lineage_tree.nwk", format_root_node=True)

if __name__ == '__main__':
    date = sys.argv[1]
    logger = logging.getLogger(__name__)
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
	                                         os.pardir))
    new_path = os.path.join(locDir, date)
    if not os.path.exists(new_path):
        os.mkdir(new_path)
        os.chdir(new_path)
    else:
        logger.warning('Outpath already exists and will be overwritten!')
    os.chdir(new_path)
    TreePath = download_tree(new_path)
    generate_files_from_tree(TreePath)
    edit_all_paths()
    df = pd.read_csv("/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/Usher_processing/"+str(date)+"/lineagePaths.txt", sep="\t")
    clade = "clade"
    file_clade = parse_tree_paths(df, clade)
    file_clade.to_csv("lineagePaths_edited_clades_"+str(date)+".tsv", sep="\t", index=False)
    file_clade.to_pickle("lineagePaths_edited_clades_"+str(date)+".pkl")
    clade = "clade_out"
    file = parse_tree_paths(df, clade)
    file.to_csv("lineagePaths_edited_"+str(date)+".tsv", sep="\t", index=False)
    file.to_pickle("lineagePaths_edited_"+str(date)+".pkl")
    create_lineage_tree_for_itol()