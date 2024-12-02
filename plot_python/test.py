import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.style.CMS)

filename = "/eos/home-j/jiehan/root/skimmed_ntuples_rui/ZGToLLG/ggfvbf_ntuples_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_2016APV201620172018.root"
data = uproot.open(filename)["tree"].arrays(library="pd")
# print(data.columns)
fig, ax = plt.subplots()
hist, bins = np.histogram(data["lly_m"], bins=100, range=(80, 180), weights=data["weightXyear"])
hep.histplot(hist, bins, histtype="fill", ax=ax)
ax.set_xlabel("Higgs Mass [GeV]")
ax.set_ylabel("Events")
plt.savefig("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/hmass_dis.png")