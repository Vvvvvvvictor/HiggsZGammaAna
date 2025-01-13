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
x_range = (100, 180, 80)
slice_var = "H_ptt"

filePath = "/eos/home-j/jiehan/root/outputs/test/zero_to_one_jet/"
bkgmc_names = ["ZGToLLG", "DYJetsToLL", "ZG2JToG2L2J"]

bkgmc_data = pd.DataFrame()
for bkgmc_name in bkgmc_names:
    fileName = filePath + bkgmc_name + ".root"
    arrays, var_not_exist = pic.read_root_file(fileName, [x_var, slice_var, "weight", "bdt_score_t", "n_jets"], "zero_to_one_jet", ["bdt_score_t>0.97"])
    bkgmc_data = bkgmc_data.append(arrays, ignore_index=True)

slice_points = bkgmc_data.query("H_ptt>-1")[slice_var].quantile(np.linspace(0, 1, 4)).values
x_edges = np.linspace(x_range[0], x_range[1], x_range[2])

for i, point in enumerate(slice_points):
    if i != len(slice_points)-1:
        bkgmc = bkgmc_data.query(f"{slice_var}>{point} & {slice_var}<{slice_points[i+1]}")
    else:
        break
    bkg_hist, _ = np.histogram(bkgmc[x_var], bins=x_edges, weights=bkgmc["weight"])
    
    fig, ax = plt.subplots()
    
    hep.histplot(bkg_hist, x_edges, label="bkgmc", histtype="fill", ax=ax)
    
    ax.set_xlabel(x_label)
    ax.set_ylabel("Events")
    ax.set_title(f"Slice in {slice_var} from {point:.2f} to {slice_points[i+1]:.2f}")
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{slice_var}_{point:.2f}-{slice_points[i+1]:.2f}.png")
    plt.close()
