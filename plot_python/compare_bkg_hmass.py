import pandas as pd
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

import os

VAR = "gamma_r9"
XLABLE = r"$R_{9}^\gamma$"
RMIN = 0.2
RMAX = 1.2
# VAR = "gamma_alliso"
# XLABLE = r"$I_{\gamma}$"
# RMIN = 0.
# RMAX = 5
BINS = 50
WEIGHT = "weight"
# PATH = "/eos/user/j/jiehan/root/outputs/"
PATH = "/eos/home-j/jiehan/data_for_norm_float_v1/"
# PATH = "/eos/home-j/jiehan/root/skimmed_ntuples/"
# "ZGToLLG" "Data"
dy1 = "ZGToLLG"
dy2 = "ZGToLLG"

def convert_root_to_hist(file, bound=None):
    type_hist = np.zeros(BINS)
    print("Reading", PATH+file+"/*.root:two_jet", "...")
    for f in os.listdir(PATH+file):
        if "root" in f:
            samples = uproot.open(PATH+file+f"/{f}:two_jet").arrays([VAR, "weight", "gamma_mvaID_WP80", "gamma_mvaID_WP90", "gamma_mvaID_WPL"], library="pd")
            samples = samples[(samples["gamma_mvaID_WP80"] == bound) & (samples["gamma_mvaID_WPL"] == 1) & (samples["gamma_mvaID_WP90"] == bound)]
            hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
            type_hist = type_hist + hist
    return type_hist, bins

# def convert_root_to_hist(file, bound=None):
#     type_hist = np.zeros(BINS)
#     print("Reading", PATH+"two_jet/"+file+".root:test", "...")
#     if "H_mass" not in VAR:
#         samples = uproot.open(PATH+"two_jet/"+file+f".root:test").arrays([VAR, "bdt_score_t", "weight", "H_mass"], library="pd")
#     else:
#         samples = uproot.open(PATH+"two_jet/"+file+".root:test").arrays([VAR, "bdt_score_t", "weight"], library="pd")   
#     if "data" not in file and "H_mass" in VAR:
#         pass
#     else:
#         samples = samples[(samples["H_mass"] < 120) | (samples["H_mass"] > 130)]
#     samples = samples[(samples["bdt_score_t"] > bound[0]) & (samples["bdt_score_t"] < bound[1])]
#     hist, bins = np.histogram(samples[VAR], bins=BINS, range=[RMIN, RMAX], weights=samples[WEIGHT])
#     type_hist = type_hist + hist

#     return type_hist, bins

hist1, bins = convert_root_to_hist(dy1, bound=1)
hist2, bins = convert_root_to_hist(dy2, bound=0)

bins = (bins[:-1] + bins[1:]) / 2
print((hist1*bins).sum()/hist1.sum())

plt.errorbar(bins, hist1/sum(hist1), yerr=np.sqrt(hist1)/sum(hist1), fmt="o", label="SR")
plt.errorbar(bins, hist2/sum(hist2), yerr=np.sqrt(hist2)/sum(hist2), fmt="o", label="CR")

# plt.yscale("log")    
plt.xlabel(XLABLE)
plt.ylabel("Events")
plt.legend()
plt.grid()
plt.tight_layout()
plt.title(f"Background component in cat.s")
plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/compare_{VAR}_Tight_Loose_ZG.png")
plt.clf()

# with open(f"{PATH}significances/bin_binaries_1D_two_jet.txt") as f:
#     lines = f.readlines()
#     bounds = [float(i) for i in lines[2].split()]
#     bounds = [0] + bounds[1:]

# for i in range(-1, len(bounds)-2):
#     bound = (bounds[i+1], bounds[i+2])
#     print(i, " ", bound)
#     hist1, bins = convert_root_to_hist(dy1, bound)
#     hist2, bins = convert_root_to_hist(dy2, bound)

#     bins = (bins[:-1] + bins[1:]) / 2
#     print((hist1*bins).sum()/hist1.sum())

#     plt.errorbar(bins, hist1/sum(hist1)+(i-1)*0.1, yerr=np.sqrt(hist1)/sum(hist1), fmt="o", label="cat{}:SM ZG".format(i))
#     plt.errorbar(bins, hist2/sum(hist2)+(i-1)*0.1, yerr=np.sqrt(hist2)/sum(hist2), fmt="o", label="cat{}:bkg".format(i))
#     # plt.errorbar(bins, hist1, yerr=np.sqrt(hist1), fmt="o", label="cat{}:bkg".format(i))
#     # plt.errorbar(bins, hist2, yerr=np.sqrt(hist2), fmt="o", label="cat{}:SM ZG".format(i))
    
#     # plt.yscale("log")    
#     plt.xlabel(XLABLE)
#     plt.ylabel("Events")
#     plt.legend()
#     plt.grid()
#     plt.tight_layout()
#     plt.title(f"Background component in cat.s")
#     plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/compare_{VAR}_ZG_bkg_{i}.png")
#     # plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/compare_{VAR}_ZG_bkg_norm.png")
#     plt.clf()