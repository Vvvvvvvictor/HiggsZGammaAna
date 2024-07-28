import uproot
import pandas as pd
import os

# Load data
data_list = ["Data", "ggH", "VBF", "ZGToLLG", "DYJetsToLL", "EWKZ2J"]
path = "/eos/user/j/jiehan/root/skimmed_ntuples/"
for data in data_list:
    file_list = os.listdir(f"{path}{data}/")
    print(data)
    for file in file_list:
        print(file)
        f = uproot.open(f"{path}{data}/{file}")
        tree = f["zero_to_one_jet"].arrays(["MET_pt"], library="pd")
        print(len(tree[tree["MET_pt"]<90]))
    print()
    