import os
import pic_template as pic
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import pandas as pd
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

x_var = "jet_1_pt"
x_label = r"$p_{T}^{j1}$ [GeV]"
x_range = (30, 130, 50)

# x_var = "jet_1_phi"
# x_label = r"$\phi^{j1}$"
# x_range = (-3.15, 3.15, 63)

y_var = "jet_1_eta"
y_label = r"$\eta^{j1}$"
y_range = (-4.7, 4.7, 47)

filePath = "/eos/home-j/jiehan/parquet/nanov12"
bkgmc_names = ["DYGto2LG_10to50", "DYGto2LG_50to100", "DYJetsToLL", "ZG2JToG2L2J"]
data_names = ["Data"]
years = ["2022preEE", "2022postEE"]

bkgmc_data = pd.DataFrame()
for bkgmc_name in bkgmc_names:
    for year in years:
        fileName = filePath + "/bkgmc/" + bkgmc_name + "_" + year + "/merged_nominal.parquet"
        arrays, var_not_exist = pic.read_parquet_file(fileName, [x_var, y_var])
        bkgmc_data = bkgmc_data.append(arrays, ignore_index=True)
        
data_data = pd.DataFrame()
for data_name in data_names:
    for year in years:
        fileName = filePath + "/data/" + data_name + "_" + year + "/merged_nominal.parquet"
        arrays, var_not_exist = pic.read_parquet_file(fileName, [x_var, y_var])
        data_data = data_data.append(arrays, ignore_index=True)
        
x_edges = np.linspace(x_range[0], x_range[1], x_range[2])
y_edges = np.linspace(y_range[0], y_range[1], y_range[2])

bkg_hist, _ , _ = np.histogram2d(bkgmc_data[x_var], bkgmc_data[y_var], bins=[x_edges, y_edges], weights=bkgmc_data["weight_central"])
data_hist, _ , _ = np.histogram2d(data_data[x_var], data_data[y_var], bins=[x_edges, y_edges], weights=data_data["weight_central"])

# set_trace()
bkg_hist = np.where(bkg_hist == 0, np.nan, bkg_hist)
ratio_hist = data_hist / bkg_hist
ratio_hist = np.nan_to_num(ratio_hist, nan=1)

fig, ax = plt.subplots()

cmap = plt.get_cmap("coolwarm")
norm = mcolor.TwoSlopeNorm(vmin=0, vcenter=1, vmax=2)

c = ax.imshow(ratio_hist.T, origin='lower', cmap=cmap, norm=norm)

cbar = fig.colorbar(c, ax=ax)
cbar.set_label("Data/MC")

ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
ax.set_title("Data/MC ratio")


ax.set_xticks((0, 10, 20, 30, 40, 50))
ax.set_xticklabels((30, 50, 70, 90, 110, 130))

# ax.set_xticks((3, 13, 23, 33, 43, 53, 63))
# ax.set_xticklabels((-3, -2, -1, 0, 1, 2, 3))

ax.set_yticks((7, 17, 27, 37, 47))
ax.set_yticklabels((-4, -2, 0, 2, 4))


plt.tight_layout()

plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/ratio_2D_{x_var}_{y_var}.png")