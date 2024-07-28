import uproot
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

# Setting
data_path = "/eos/home-j/jiehan/data_for_norm_float_v1/output/two_jet/"
region = "two_jet"
variable = "H_mass"
xlabel = r"$m_{ll\gamma}$ [GeV/c$^2$]"
xmin, xmax, bins = 100, 180, 80
data_list = ["Data"]
bkg_list = ["ZG2JToG2L2J", "ZGToLLG"]

selection = "bdt_score_t > 0.92"

# Load data
data_hists = []
for data in data_list:
    print("reading: ", data_path + data + "_fake.root")
    data_hist = np.histogram(uproot.open(data_path + data + "_fake.root")[region].arrays(library="pd").query(selection)[variable], bins=bins, range=(xmin, xmax))[0]
    data_hists.append(data_hist)
        
bkg_hists = []
for bkg in bkg_list:
    print("reading: ", data_path + bkg + "_fake.root")
    bkg_hist = np.histogram(uproot.open(data_path + bkg + "_fake.root")[region].arrays(library="pd").query(selection)[variable], bins=bins, range=(xmin, xmax))[0]
    bkg_hists.append(bkg_hist)
    
# Plot
plt.figure()
hep.histplot(bkg_hists, bins=np.linspace(xmin, xmax, bins+1), label=bkg_list, stack=True, histtype="fill")
plt.errorbar(np.linspace(xmin, xmax, bins), data_hists[0], fmt='o', label="Data", yerr=np.sqrt(data_hists[0]), color='k')

plt.xlabel(xlabel)
plt.ylabel("Events")
plt.legend(fontsize=16, loc="best")
plt.title(f"Control region: {xmin} < {variable} < {xmax}")
plt.xlim(xmin, xmax)
plt.grid()
plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/compare_control_region_{variable}.png")