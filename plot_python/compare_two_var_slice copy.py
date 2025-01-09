import pandas as pd
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

import os
import pic_template as pic

x_var = "jet_1_eta"
x_label = r"$\eta^{j1}$"
x_range = (-4.7, 4.7, 47)
slice_var = "jet_1_pt"

filePath = "/eos/home-j/jiehan/parquet/nanov12"
bkgmc_names = ["DYGto2LG_10to50", "DYGto2LG_50to100", "DYJetsToLL", "ZG2JToG2L2J"]
data_names = ["Data"]
years = ["2022preEE", "2022postEE"]

bkgmc_data = pd.DataFrame()
for bkgmc_name in bkgmc_names:
    for year in years:
        fileName = filePath + "/bkgmc/" + bkgmc_name + "_" + year + "/merged_nominal.parquet"
        arrays, var_not_exist = pic.read_parquet_file(fileName, [x_var, slice_var])
        bkgmc_data = bkgmc_data.append(arrays, ignore_index=True)
        
data_data = pd.DataFrame()
for data_name in data_names:
    for year in years:
        fileName = filePath + "/data/" + data_name + "_" + year + "/merged_nominal.parquet"
        arrays, var_not_exist = pic.read_parquet_file(fileName, [x_var, slice_var])
        data_data = data_data.append(arrays, ignore_index=True)
        
slice_points = data_data.query("jet_1_pt>0")[slice_var].quantile(np.linspace(0, 1, 10)).values
x_edges = np.linspace(x_range[0], x_range[1], x_range[2])
    
# fig, ax = plt.subplots()

for i, point in enumerate(slice_points):
    if i != len(slice_points)-1:
        bkgmc = bkgmc_data.query(f"{slice_var}>{point} & {slice_var}<{slice_points[i+1]}")
        data = data_data.query(f"{slice_var}>{point} & {slice_var}<{slice_points[i+1]}")
    else:
        break
    bkg_hist, _ = np.histogram(bkgmc[x_var], bins=x_edges, weights=bkgmc["weight_central"])
    data_hist, _ = np.histogram(data[x_var], bins=x_edges, weights=data["weight_central"])
    
    bkg_hist = np.where(bkg_hist == 0, np.nan, bkg_hist)
    ratio_hist = data_hist / bkg_hist
    ratio_hist = np.nan_to_num(ratio_hist, nan=1)
    
    fig, ax = plt.subplots()
    
    ax.errorbar(x_edges[:-1], ratio_hist, yerr=np.sqrt(data_hist)/bkg_hist, fmt='o', label=f"{point:.2f}-{slice_points[i+1]:.2f}", markersize=10)
    
    ax.set_xlabel(x_label)
    ax.set_ylabel("Data/MC")
    ax.set_title(f"Slice in {slice_var} from {point:.2f} to {slice_points[i+1]:.2f}")
    ax.legend()

    ax.set_ylim(0, 2)
        
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{slice_var}_{x_var}_{point:.2f}-{slice_points[i+1]:.2f}.png")
    plt.clf()