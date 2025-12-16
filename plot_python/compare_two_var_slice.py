import pandas as pd
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
plt.style.use(hep.style.CMS)

import os
import pic_template as pic

x_var = "H_mass"
x_label = r"$m_{\ell\ell\gamma}$[GeV]"
x_range = (115, 180, 65)
slice_var = "bdt_score_t"

filePath = "/eos/home-j/jiehan/root/outputs/test/zero_to_one_jet/"
bkgmc_names = ["ZGToLLG", "DYJetsToLL_ext", "ZG2JToG2L2J"]

bkgmc_data = pd.DataFrame()
for bkgmc_name in bkgmc_names:
    fileName = filePath + bkgmc_name + ".root"
    arrays, var_not_exist = pic.read_root_file(fileName, [x_var, slice_var, "weight", "bdt_score_t", "n_jets"], "zero_to_one_jet", [])
    bkgmc_data = bkgmc_data.append(arrays, ignore_index=True)

slice_points = bkgmc_data.query("H_mass>115")[slice_var].quantile(np.linspace(0, 1, 10)).values
x_edges = np.linspace(x_range[0], x_range[1], x_range[2])

fig, ax = plt.subplots()

for i, point in enumerate(slice_points):
    if i != len(slice_points)-1:
        bkgmc = bkgmc_data.query(f"{slice_var}>{point} & {slice_var}<{slice_points[i+1]}")
    else:
        break
    bkg_hist, _ = np.histogram(bkgmc[x_var], bins=x_edges, weights=bkgmc["weight"])
    
    hep.histplot(bkg_hist, x_edges, label=f"{point:.3f} < {slice_var} < {slice_points[i+1]:.3f}", histtype="step", ax=ax)

ax.set_xlabel(x_label)
ax.set_ylabel("Events")
ax.set_title(f"Slices in {slice_var}")
ax.legend()
plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{slice_var}_all_slices.png")
plt.close()
