import pandas as pd
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
import os


VAR = "Z_mass"
XLABLE = r"$M_{Z}$"
BINS = 500
RMIN = 75
RMAX = 105
WEIGHT = "weight_central"
PATH = "/eos/user/j/jiehan/parquet/nanov9/signal/"

mc_file_list = [
    "ggH",
    "VBFH",
    "WplusH",
    "WminusH",
    "ZH",
    "ttH"
]

mc_hist_list, data_hist_list = [], []

def convert_parquet_to_hist(file_list, hist_list):
    for data in file_list:
        hist = np.zeros(BINS)
        for sub_dir in os.listdir(PATH):
            if data in sub_dir:
                print("Processing", PATH+sub_dir+"/merged_nominal.parquet", "...")
                samples = pd.read_parquet(PATH+sub_dir+"/merged_nominal.parquet")
                # set_trace()
                # samples = samples[samples["n_jets"]>1]
                # if "weight_central_no_lumi" in samples.keys():
                #     print(data, ":", samples["weight_central"][0]/samples["weight_central_no_lumi"][0], sep="\n")
                # if "DY" in data:
                #     # print(sum(samples[samples["n_iso_photons"]>0][WEIGHT]))
                #     # print(sum(samples[samples["n_iso_photons"]==0][WEIGHT]))
                #     samples = samples[samples["n_iso_photons"]==0]
                #     # samples[WEIGHT]*=3
                # if "ZG" in data:
                #     # print(sum(samples[samples["n_iso_photons"]>0][WEIGHT]))
                #     # print(sum(samples[samples["n_iso_photons"]==0][WEIGHT]))
                #     samples = samples[samples["n_iso_photons"]>0]
                #     # samples[WEIGHT]*=3
                # if "H_mass" in VAR and "data" in data:
                #     samples = samples[(samples["H_mass"]<122) | (samples["H_mass"]>128)]
                sub_hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
                hist = hist + sub_hist
        hist_list.append(hist)

def convert_root_to_hist(file_list, hist_list):
    for type in file_list:
        type_hist = np.zeros(BINS)
        for data in type:
            print("Reading", PATH+data+"/two_jet_run2.root:two_jet", "...")
            samples = uproot.open(PATH+data+"/two_jet_run2.root:two_jet").arrays(library="pd")
            if ("H_mass" in VAR and "data" in data) or "H_mass" not in VAR:
                samples = samples[(samples["H_mass"]<122) | (samples["H_mass"]>128)]
            print(samples["weight"].sum())
            hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
            type_hist = type_hist + hist
        hist_list.append(type_hist)

convert_parquet_to_hist(mc_file_list, mc_hist_list)

hep.style.use("CMS")
fig, ax = plt.subplots(figsize=(10, 10), dpi=400)

bins = np.linspace(RMIN, RMAX, BINS+1)

hep.histplot(
    mc_hist_list, 
    bins=bins, 
    histtype='fill', 
    color=['b','g', 'm', 'y', 'r', 'c'],
    ax=ax,
    label=["ggH", "VBF", r"$W^{-}$H", r"$W^{+}$H", "ZH", r"$t\bar{t}$H"],
    stack=True
)

plt.plot(91.1876, 0, 91.1876, 2, linestyle='--', color='black', linewidth=3)

ax.set_xlabel("{} (GeV)".format(XLABLE), fontsize=25)
ax.set_ylabel("Events / {:.3f} GeV".format((RMAX-RMIN)/BINS), fontsize=25)
ax.set_xlim(RMIN, RMAX)
ax.legend()

plt.savefig("pic/2jet_{}.png".format(VAR))
