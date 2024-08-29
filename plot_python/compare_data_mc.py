import pandas as pd
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace


# VAR = "H_mass"
# XLABLE = r"$M_{ll\gamma}$"
# BINS = 40
# RMIN = 100
# RMAX = 180
VAR = "n_jets"
XLABLE = r"$n_{j}$"
BINS = 5
RMIN = 1
RMAX = 6
WEIGHT = "weight"
PATH = "/eos/user/j/jiehan/root/skimmed_ntuples/"

bkg_mc_file_list = [
    ["ZGToLLG"],
    ["DYJetsToLL"],
    ["EWKZ2J"],
    ["ZG2JToG2L2J"],
    ["TT"],
    ["TTGJets"],
    ["ttWJets", "ttZJets"],
    ["WW", "WZ", "ZZ"]
]
data_file_list = [
    ["Data"]
]

bkg_mc_hist_list, data_hist_list = [], []

def convert_parquet_to_hist(file_list, hist_list):
    for data in file_list:
        samples = pd.read_parquet(data)
        # samples = samples[samples["n_jets"]>1]
        # if "weight_central_no_lumi" in samples.keys():
        #     print(data, ":", samples["weight_central"][0]/samples["weight_central_no_lumi"][0], sep="\n")
        if "DY" in data:
            # print(sum(samples[samples["n_iso_photons"]>0][WEIGHT]))
            # print(sum(samples[samples["n_iso_photons"]==0][WEIGHT]))
            samples = samples[samples["n_iso_photons"]==0]
            # samples[WEIGHT]*=3
        # if "ZG" in data:
        #     # print(sum(samples[samples["n_iso_photons"]>0][WEIGHT]))
        #     # print(sum(samples[samples["n_iso_photons"]==0][WEIGHT]))
        #     samples = samples[samples["n_iso_photons"]>0]
        #     # samples[WEIGHT]*=3
        if "H_mass" in VAR and "data" in data:
            samples = samples[(samples["H_mass"]<122) | (samples["H_mass"]>128)]
        hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
        hist_list.append(hist)

def convert_root_to_hist(file_list, hist_list):
    for type in file_list:
        type_hist = np.zeros(BINS)
        for data in type:
            for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
                print("Reading", PATH+data+f"/{year}.root:inclusive", "...")
                samples = uproot.open(PATH+data+f"/{year}.root:inclusive").arrays([VAR, "H_mass", "weight"] if "H_mass" not in VAR else [VAR, "weight"], library="pd")
                if ("H_mass" in VAR and "Data" in data) or "H_mass" not in VAR:
                    samples = samples[(samples["H_mass"]<120) | (samples["H_mass"]>130)]
                print(samples["weight"].sum())
                hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
                type_hist = type_hist + hist
        hist_list.append(type_hist)

# convert_parquet_to_hist(data_list, data_hist_list)
# convert_parquet_to_hist(bkg_mc_list, bkg_mc_hist_list)

convert_root_to_hist(data_file_list, data_hist_list)
convert_root_to_hist(bkg_mc_file_list, bkg_mc_hist_list)

hep.style.use("CMS")
fig, ax = plt.subplots(figsize=(10, 10), dpi=400)

bins = np.linspace(RMIN, RMAX, BINS+1)
hep.histplot(
    data_hist_list, 
    bins=bins, 
    stack=True, 
    histtype='errorbar', 
    color='r',
    ax=ax,
    label=["Data"]
)
hep.histplot(
    bkg_mc_hist_list, 
    bins=bins, 
    histtype='fill', 
    color=['b','g'],
    ax=ax,
    label=["SM ZG", "DYJets", "EWK Z+Jets", "EWK ZG", "TT", "TTG+Jets", "TTVJets", "Diboson"],
    stack=True
)

ax.set_xlabel("{} (GeV)".format(XLABLE), fontsize=25)
ax.set_ylabel("Events / {:.0f} GeV".format((RMAX-RMIN)/BINS), fontsize=25)
ax.set_xlim(RMIN, RMAX)
ax.legend()

plt.savefig("pic/data_vs_mc_2jet_{}.png".format(VAR))
