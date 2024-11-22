import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.style.CMS)

bkg_names = ["ZGToLLG", "DYJetsToLL", "ZG2JToG2L2J"] #, "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
bkg_hist = []
colors = [mcolors.CSS4_COLORS["red"], mcolors.CSS4_COLORS["green"], mcolors.CSS4_COLORS["blue"]] #, mcolors.CSS4_COLORS["orange"], mcolors.CSS4_COLORS["purple"], mcolors.CSS4_COLORS["brown"], mcolors.CSS4_COLORS["pink"], mcolors.CSS4_COLORS["gray"], mcolors.CSS4_COLORS["olive"], mcolors.CSS4_COLORS["cyan"], mcolors.CSS4_COLORS["lime"], mcolors.CSS4_COLORS["teal"], mcolors.CSS4_COLORS["magenta"]]

with open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/best_boundaries.txt", "r") as f:
    n, m = list(map(float, f.readline().split(",")))
for bkg in bkg_names:
    cat_data = []
    data = uproot.open(f"/eos/user/j/jiehan/root/outputs/two_jet/{bkg}.root")["two_jet"].arrays(["lly_m", "bdt_score_t", "vbf_score_t", "w_lumiXyear"], library="pd")
    cat_data.append(data.query(f"bdt_score_t<{n} & vbf_score_t<{m}"))
    cat_data.append(data.query(f"bdt_score_t<{n} & vbf_score_t>{m}"))
    cat_data.append(data.query(f"bdt_score_t>{n} & vbf_score_t<{m}"))
    cat_data.append(data.query(f"bdt_score_t>{n} & vbf_score_t>{m}"))
    for data_c in cat_data:
        hist, bins = np.histogram(data_c["lly_m"], bins=100, weights=data_c["w_lumiXyear"], range=(100, 180))
        bkg_hist.append(hist)

for i in range(4):
    fig, ax = plt.subplots()
    plot_hist = [bkg_hist[i] for i in range(i, len(bkg_hist), 4)]
    hep.histplot(plot_hist, bins, stack=True, label=bkg_names, histtype="fill", ax=ax, color=colors)
    ax.set_xlim(100, 180)
    ax.set_xlabel("Higgs Mass [GeV]")
    ax.set_ylabel("Events")
    ax.legend()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/2D_hmass_dis_{i}.png")
    
fig, ax = plt.subplots()
plot_hist = [np.sum(bkg_hist[i:i+4], axis=0) for i in range(0, len(bkg_hist), 4)]
hep.histplot(plot_hist, bins, stack=True, label=bkg_names, histtype="fill", ax=ax, color=colors)
ax.set_xlim(100, 180)
ax.set_xlabel("Higgs Mass [GeV]")
ax.set_ylabel("Events")
ax.legend()
plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/2D_hmass_dis_all.png")