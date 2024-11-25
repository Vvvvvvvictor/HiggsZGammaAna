import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.style.CMS)

bkg_names = ["ZGToLLG", "DYJetsToLL", "ZG2JToG2L2J"]
signal_names = ["ggH_M125", "VBF_M125"]
data_names = ["Data"]
bkg_hist, sig_hist, data_hist = [], [], []
colors = [mcolors.CSS4_COLORS["green"], mcolors.CSS4_COLORS["orange"], mcolors.CSS4_COLORS["blue"]] #, mcolors.CSS4_COLORS["orange"], mcolors.CSS4_COLORS["purple"], mcolors.CSS4_COLORS["brown"], mcolors.CSS4_COLORS["pink"], mcolors.CSS4_COLORS["gray"], mcolors.CSS4_COLORS["olive"], mcolors.CSS4_COLORS["cyan"], mcolors.CSS4_COLORS["lime"], mcolors.CSS4_COLORS["teal"], mcolors.CSS4_COLORS["magenta"]]

with open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/best_boundaries.txt", "r") as f:
    p = list(map(float, f.readline().strip().split(",")))
    significances = list(map(float, f.readline().split(",")))
bin_edges_x = [0, 1]
bin_edges_y = [0, 1]
for i in p[:-1]:
    bin_edges_x.insert(-1, i)
bin_edges_y.insert(-1, p[-1])
for bkg in bkg_names:
    cat_data = []
    data = uproot.open(f"/eos/user/j/jiehan/root/outputs/two_jet/{bkg}.root")["two_jet"].arrays(["H_mass", "bdt_score_t", "vbf_score_t", "weight"], library="pd")
    for i in range(len(bin_edges_x)-1):
        for j in range(len(bin_edges_y)-1):
            data_c = data.query(f"bdt_score_t>{bin_edges_x[i]} & bdt_score_t<{bin_edges_x[i+1]} & vbf_score_t>{bin_edges_y[j]} & vbf_score_t<{bin_edges_y[j+1]}")
            cat_data.append(data.query(f"bdt_score_t>{bin_edges_x[i]} & bdt_score_t<{bin_edges_x[i+1]} & vbf_score_t>{bin_edges_y[j]} & vbf_score_t<{bin_edges_y[j+1]}"))
    for data_c in cat_data:
        hist, bins = np.histogram(data_c["H_mass"], bins=80, weights=data_c["weight"], range=(100, 180))
        bkg_hist.append(hist)
        
for sig in signal_names:
    cat_data = []
    data = uproot.open(f"/eos/user/j/jiehan/root/outputs/two_jet/{sig}.root")["two_jet"].arrays(["H_mass", "bdt_score_t", "vbf_score_t", "weight"], library="pd")
    for i in range(len(bin_edges_x)-1):
        for j in range(len(bin_edges_y)-1):
            data_c = data.query(f"bdt_score_t>{bin_edges_x[i]} & bdt_score_t<{bin_edges_x[i+1]} & vbf_score_t>{bin_edges_y[j]} & vbf_score_t<{bin_edges_y[j+1]}")
            cat_data.append(data.query(f"bdt_score_t>{bin_edges_x[i]} & bdt_score_t<{bin_edges_x[i+1]} & vbf_score_t>{bin_edges_y[j]} & vbf_score_t<{bin_edges_y[j+1]}"))
    for data_c in cat_data:
        hist, bins = np.histogram(data_c["H_mass"], bins=80, weights=data_c["weight"], range=(100, 180))
        sig_hist.append(hist)
        
for data in data_names:
    cat_data = []
    data = uproot.open(f"/eos/user/j/jiehan/root/outputs/two_jet/{data}.root")["two_jet"].arrays(["H_mass", "bdt_score_t", "vbf_score_t", "weight"], library="pd").query("H_mass < 120 | H_mass > 130")
    for i in range(len(bin_edges_x)-1):
        for j in range(len(bin_edges_y)-1):
            data_c = data.query(f"bdt_score_t>{bin_edges_x[i]} & bdt_score_t<{bin_edges_x[i+1]} & vbf_score_t>{bin_edges_y[j]} & vbf_score_t<{bin_edges_y[j+1]}")
            cat_data.append(data.query(f"bdt_score_t>{bin_edges_x[i]} & bdt_score_t<{bin_edges_x[i+1]} & vbf_score_t>{bin_edges_y[j]} & vbf_score_t<{bin_edges_y[j+1]}"))
    for data_c in cat_data:
        hist, bins = np.histogram(data_c["H_mass"], bins=80, weights=data_c["weight"], range=(100, 180))
        data_hist.append(hist)

num_cat = (len(bin_edges_x) - 1) * (len(bin_edges_y) - 1)
cat_num = []
for i in range(len(bin_edges_x)-1):
    for j in range(len(bin_edges_y)-1):
        cat_num.append(j*(len(bin_edges_x)-1) + i)
for i in range(num_cat):
    fig, ax = plt.subplots()
    plot_hist = [bkg_hist[i] for i in range(i, len(bkg_hist), num_cat)]
    plot_sig = np.sum([sig_hist[i] for i in range(i, len(sig_hist), num_cat)], axis=0)
    plot_data = np.sum([data_hist[i] for i in range(i, len(data_hist), num_cat)], axis=0)
    hep.histplot(plot_hist, bins, stack=True, label=bkg_names, histtype="fill", ax=ax, color=colors)
    ax.errorbar((bins[:-1] + bins[1:]) / 2, plot_data, yerr=np.sqrt(plot_data), fmt="o", label="data", color="black")
    ax.plot((bins[:-1] + bins[1:]) / 2, plot_sig, label="signal", color="red")
    ax.set_xlim(100, 180)
    ax.set_xlabel("Higgs Mass [GeV]")
    ax.set_ylabel("Events")
    ax.annotate(f"Significance: {significances[i]:.4f}", xy=(0.1, 0.95), xycoords="axes fraction")
    ax.legend()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/2D_hmass_dis_{cat_num[i]}.png")
    
fig, ax = plt.subplots()
plot_hist = [np.sum(bkg_hist[i:i+num_cat], axis=0) for i in range(0, len(bkg_hist), num_cat)]
hep.histplot(plot_hist, bins, stack=True, label=bkg_names, histtype="fill", ax=ax, color=colors)
plot_sig = np.sum(sig_hist, axis=0)
ax.errorbar((bins[:-1] + bins[1:]) / 2, plot_data, yerr=np.sqrt(plot_data), fmt="o", label="data", color="black")
ax.plot((bins[:-1] + bins[1:]) / 2, plot_sig, label="signal", color="red")
ax.set_xlim(100, 180)
ax.set_xlabel("Higgs Mass [GeV]")
ax.set_ylabel("Events")
ax.legend()
plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/2D_hmass_dis_all.png")