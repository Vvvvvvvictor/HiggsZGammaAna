import uproot
import pandas as pd
import numpy as np
from pdb import set_trace

path = "/eos/home-j/jiehan/root/outputs/two_jet"
# path = "/eos/home-j/jiehan/data_for_norm_float_v1/output_WP80/zero_to_one_jet/"
filename = "data_driven_bkg_v3"
tree = "two_jet"
# filename = "DYJets_01J"
# tree = "zero_to_one_jet"
# filenames = ["DYJets_01J", "ZGToLLG", "ZG2JToG2L2J"]
# trees = ["zero_to_one_jet", "zero_to_one_jet", "zero_to_one_jet"]

# boundaries = [0.2838, 0.4558, 0.5797, 0.7070, 1.00]
boundaries = [float(i) for i in open("/eos/home-j/jiehan/root/outputs/significances/bin_binaries_1D_two_jet.txt", "r").readlines()[2].split()[1:]]

# data = pd.DataFrame()
# for i, filename in enumerate(filenames):
#     print(f"path: {path}/{filename}:{trees[i]}")
#     data = uproot.open(f"{path}/{filename}.root:{trees[i]}").arrays(library="pd")
#     for i in range(len(boundaries)-1):
#         temp = data.query(f"bdt_score >= {boundaries[i]} and bdt_score < {boundaries[i+1]}")
#         error_r = (sum(temp["weight"]))**0.5 / (sum((temp["weight"].to_numpy())**2))**0.5
#         weighted = temp['weight'].sum()
#         # if "data_driven" in filename:
#         #     temp = temp.query("weight > 0.01")
#         unweighted = temp.shape[0]
#         print("minus: ", temp[temp["weight"]<0].shape[0])
#         for num in np.unique(temp["weight"]):
#             print(f"num: {num}, shape: {temp[temp['weight'] == num].shape[0]}")
#         print(f"bin {i+1}: unweighted: {unweighted}, weighted: {weighted}, error_r: {error_r}")

print(f"path: {path}/{filename}:{tree}")
data = uproot.open(f"{path}/{filename}.root:{tree}").arrays(library="pd")
output_file = open(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/output_{filename}.txt", "w")
for i in range(len(boundaries)-1):
    temp = data.query(f"bdt_score_t >= {boundaries[i]} and bdt_score_t < {boundaries[i+1]}")
    error_r = (sum(temp["weight_cor"]))**0.5 / (sum((temp["weight_cor"].to_numpy())**2))**0.5
    weighted = temp['weight_cor'].sum()
    # if "data_driven" in filename:
    #     temp = temp.query("weight > 0.01")
    unweighted = temp.shape[0]
    # print("minus: ", temp[temp["weight_cor"]<0].shape[0])
    # for num in np.unique(temp["weight_cor"]):
    #     print(f"num: {num}, shape: {temp[temp['weight_cor'] == num].shape[0]}")
    print(f"bin {i+1}: unweighted: {unweighted}, weighted: {weighted}, error_r: {error_r}")
    output_file.write(f"bin {i+1}: unweighted: {unweighted}, weighted: {weighted}, error_r: {error_r}\n")