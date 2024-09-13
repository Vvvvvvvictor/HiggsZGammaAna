import pandas as pd
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

import os
import pic_template as pic

# VAR = "gamma_eta"
# XLABLE = r"$\eta_\gamma$"
# RMIN = -2.5
# RMAX = 2.5
VAR = "H_ptt"
XLABLE = r"$m_{\ell\ell\gamma}$"
RMIN = 0
RMAX = 50
slice_var = "H_mass"
BINS = 50
WEIGHT = "weight"
# PATH = "/eos/user/j/jiehan/root/outputs/two_jet"
# PATH = "/eos/home-j/jiehan/data_for_norm_float_v1/"
PATH = "/eos/home-j/jiehan/root/skimmed_ntuples/"
# "ZGToLLG" "Data"
TREE = "two_jet"
sample = "Data"
selection = []

def convert_root_to_hist(sample, boundaries):
    data = pd.DataFrame()
    for file in os.listdir(PATH+sample):
        if file.endswith(".root"): 
            data_temp, _ = pic.read_file(PATH+sample+"/"+file, var=[VAR, "weight", slice_var], tree=TREE, selections=selection)
            data = pd.concat([data, data_temp])
    hist_list, hist_err_list = [], []
    for i in range(len(boundaries)-1):
        bounds = boundaries[i:i+2]
        data_slice = data[(data[slice_var] > bounds[0]) & (data[slice_var] < bounds[1])]
        hist, bins = np.histogram(data_slice[VAR], bins=BINS, range=[RMIN, RMAX], weights=data_slice[WEIGHT])
        hist_err, _ = np.histogram(data_slice[VAR], bins=BINS, range=[RMIN, RMAX], weights=data_slice[WEIGHT]**2)
        hist_list.append(hist)
        hist_err_list.append(hist_err)
    return hist_list, bins, hist_err_list

data = pd.DataFrame()
for file in os.listdir(PATH+sample):
    if file.endswith(".root"): 
        data_temp, _ = pic.read_file(PATH+sample+"/"+file, var=[VAR, "weight", slice_var], tree=TREE, selections=selection)
        data = pd.concat([data, data_temp])
# get the boundaries that are used to slice the data, devide samples according to the percentage of slice var
boundaries = data[slice_var].quantile(np.linspace(0, 1, 6)).values
hist_list, bins, hist_err_list = convert_root_to_hist(sample, boundaries)

bins = (bins[:-1] + bins[1:]) / 2

for i, hist in enumerate(hist_list):
    plt.errorbar(bins, hist, yerr=np.sqrt(hist_err_list[i]), fmt="o", label="{:.2f} to {:.2f}".format(boundaries[i], boundaries[i+1]))
    # plt.plot(bins, hist)

# plt.yscale("log")    
plt.xlabel(XLABLE)
plt.ylabel("Events")
plt.legend()
plt.grid()
plt.tight_layout()
plt.title(f"{VAR} shape in {slice_var} percentage regions")
plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/compare_data_{VAR}_{slice_var}.png")
plt.clf()
