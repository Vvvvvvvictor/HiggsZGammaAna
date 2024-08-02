import uproot
import pandas as pd
import numpy as np
from pdb import set_trace

path = "/eos/home-j/jiehan/root/outputs/zero_to_one_jet"
# filename = "data_driven_bkg"
# tree = "test"
# filename = "DYJets_01J"
# tree = "zero_to_one_jet"
filenames = ["DYJets_01J", "ZGToLLG", "ZG2JToG2L2J"]
trees = ["zero_to_one_jet", "zero_to_one_jet", "zero_to_one_jet"]

boundaries = [0.2838, 0.4558, 0.5797, 0.7070, 1.00]

data = pd.DataFrame()
for i, filename in enumerate(filenames):
    print(f"path: {path}/{filename}:{trees[i]}")
    data = uproot.open(f"{path}/{filename}.root:{trees[i]}").arrays(library="pd")
    for i in range(len(boundaries)-1):
        temp = data.query(f"bdt_score >= {boundaries[i]} and bdt_score < {boundaries[i+1]}")
        error_r = (sum(temp["weight"]))**0.5 / (sum((temp["weight"].to_numpy())**2))**0.5
        weighted = temp['weight'].sum()
        # if "data_driven" in filename:
        #     temp = temp.query("weight > 0.01")
        unweighted = temp.shape[0]
        print("minus: ", temp[temp["weight"]<0].shape[0])
        for num in np.unique(temp["weight"]):
            print(f"num: {num}, shape: {temp[temp['weight'] == num].shape[0]}")
        print(f"bin {i+1}: unweighted: {unweighted}, weighted: {weighted}, error_r: {error_r}")