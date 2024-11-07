import awkward as ak
import pandas as pd
import uproot
import matplotlib.pyplot as plt
from pdb import set_trace

# data = pd.read_parquet("/eos/home-j/jiehan/parquet/nanov9/cutflow/ggH_M125_2017/merged_nominal.parquet")
# set_trace()
# log = pd.read_csv("jets.txt", sep=" ", header=None)

# log_data = pd.DataFrame()
# for i in range(0, len(log)):
#     temp = data[(data["run"] == log[0][i]) & (data["luminosityBlock"] == log[1][i]) & (data["event"] == log[2][i])]
#     if temp.empty:
#         print(log.loc[:,i])
#     else:
#         log_data = pd.concat([log_data, temp])
        
# index = log_data.index.to_list()
# data = data.drop(index)
# exo_data = data[data["n_jets"]>0]
# set_trace()

# data = uproot.open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/hzgml/ggH_M125/2017.root")["inclusive"].arrays(["lep_cos_theta", "Z_cos_theta", "lep_phi", "z_mumu", "z_ee", "weight"],library="pd")
# cut1 = data["lep_cos_theta"] > 0.5
# cut2 = data["Z_cos_theta"] > 0.5
# cut3 = data["lep_phi"] > 0.5

# mu = data["z_mumu"] == 1
# ee = data["z_ee"] == 1

# print("ee:")
# print(sum(data[ee & cut1]["weight"]))
# print(sum(data[ee & cut2]["weight"]))
# print(sum(data[ee & cut3]["weight"]))
# print(sum(data[ee & cut1 & cut2 & cut3]["weight"]))
# print("mu:")
# print(sum(data[mu & cut1]["weight"]))
# print(sum(data[mu & cut2]["weight"]))
# print(sum(data[mu & cut3]["weight"]))
# print(sum(data[mu & cut1 & cut2 & cut3]["weight"]))

data = uproot.open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA/13481923-E09A-FE40-97E8-C04572175215.root:Events").arrays(["Photon_electronVeto", "Photon_pt", "Photon_isScEtaEB", "Photon_isScEtaEE"], library="ak")
# data = uproot.open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/HiggsDNA/0053002A-3866-0F4A-912C-12E0440839CEskimmed.root:Events").arrays(["Photon_electronVeto", "Photon_pt", "Photon_isScEtaEB", "Photon_isScEtaEE"], library="ak")
print(ak.sum(data["Photon_electronVeto"])/len(ak.flatten(data["Photon_electronVeto"])))
# fig = plt.figure()
# plt.hist(data["Photon_electronVeto"].flatten(), bins=100, range=(0,1))
# plt.savefig("test.png")