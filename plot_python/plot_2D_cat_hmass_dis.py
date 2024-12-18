import uproot
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.style.CMS)

from pdb import set_trace

data = uproot.open("output_with_categories.root")["test"].arrays(library="pd")
bkg = data.query("label == 0")
sig = data.query("label == 1")
if os.path.exists("figs") == False:
    os.mkdir("figs")
for i in range(max(data["category"])+1):
    bkg_i = bkg.query("category == {}".format(i))
    sig_i = sig.query("category == {}".format(i))
    bkg_hist, bins = np.histogram(bkg_i["H_mass"], bins=80, range=(100, 180), weights=bkg_i["weight"])
    bkg_error, _ = np.histogram(bkg_i["H_mass"], bins=80, range=(100, 180), weights=bkg_i["weight"]**2)
    sig_hist, _ = np.histogram(sig_i["H_mass"], bins=80, range=(100, 180), weights=sig_i["weight"])
    pos = (bins[:-1] + bins[1:]) / 2
    sig_err, _ = np.histogram(sig_i["H_mass"], bins=80, range=(100, 180), weights=sig_i["weight"]**2)
    
    fig, ax = plt.subplots()
    hep.histplot(bkg_hist, bins, yerr=np.sqrt(bkg_error), label="bkg", histtype="fill", color="green")
    ax.errorbar(pos, sig_hist, yerr=np.sqrt(sig_err), fmt="o", label="sig", color="red")
    ax.set_xlabel(r"$m_{ll\gamma}$ [GeV]")
    ax.set_ylabel("Events")
    ax.legend()
    plt.savefig("figs/output_{}.png".format(i))

socre_boundaries = np.percentile(sig["score"], np.linspace(0, 100, 101))
def get_significance(sig, bkg, score):
    if bkg > 0:
        return np.sqrt(2*((sig+bkg)*np.log(1+sig/bkg)-sig))
significances = []
for i in range(len(socre_boundaries)-1):
    sig_i = sig.query("score > {} and score <= {} and H_mass > 120 and H_mass < 130".format(socre_boundaries[i], socre_boundaries[i+1]))
    bkg_i = bkg.query("score > {} and score <= {} and H_mass > 120 and H_mass < 130".format(socre_boundaries[i], socre_boundaries[i+1]))
    significances.append(get_significance(sig_i["weight"].sum(), bkg_i["weight"].sum(), bkg_i["weight"].sum()))
fig, ax = plt.subplots()
ax.scatter(socre_boundaries[:-1], significances, s=25)
ax.set_xlabel("score")
ax.set_ylabel("significance")
plt.savefig("figs/significance.png")

dirpath = '/eos/home-j/jiehan/root/outputs/two_jet_2D_categories'
bkg_names = ["ZGToLLG", "DYJetsToLL", "ZG2JToG2L2J"]
sig_names = ["ggH_M125", "VBF_M125", "ZH_M125", "WplusH_M125", "WminusH_M125", "ZH_M125", "ttH_M125"]
data_names = ["Data"]

bkg_data, sig_data = [], []
for name in bkg_names:
    bkg_data.append(uproot.open(f"{dirpath}/{name}.root")["two_jet"].arrays(library="pd"))
for name in sig_names:
    sig_data.append(uproot.open(f"{dirpath}/{name}.root")["two_jet"].arrays(library="pd"))
data_data = uproot.open(f"{dirpath}/Data.root")["two_jet"].arrays(library="pd").query("H_mass < 120 | H_mass > 130")

with open("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/models/best_thresholds.txt", "r") as f:
    boundaries = list(map(float, f.readline().strip().split(" ")))
    print(boundaries)
    
for i in range(len(boundaries)-1):
    selections = f"bdt_score_t>{boundaries[i]} & bdt_score_t<{boundaries[i+1]} & H_mass>100 & H_mass<180"
    bkg_hists = [np.histogram(data.query(selections)["H_mass"], bins=80, range=(100, 180), weights=data.query(selections)["weight"])[0] for data in bkg_data]
    bkg_err_hists = [np.histogram(data.query(selections)["H_mass"], bins=80, range=(100, 180), weights=data.query(selections)["weight"]**2)[0] for data in bkg_data]
    sig_hist = np.histogram(np.concatenate([data.query(selections)["H_mass"] for data in sig_data]), bins=80, range=(100, 180), weights=np.concatenate([data.query(selections)["weight"] for data in sig_data]))[0]
    VBF_hist = np.histogram(sig_data[1].query(selections)["H_mass"], bins=80, range=(100, 180), weights=sig_data[1].query(selections)["weight"])[0]
    data_hist = np.histogram(data_data.query(selections)["H_mass"], bins=80, range=(100, 180), weights=data_data.query(selections)["weight"])[0]
    
    fig, ax = plt.subplots()
    bins = np.linspace(100, 180, 81)
    pos = (bins[:-1] + bins[1:]) / 2
    
    hep.histplot(bkg_hists, bins, yerr=np.sqrt(bkg_err_hists), label=bkg_names, histtype="fill", stack=True, ax=ax)
    ax.errorbar(pos, data_hist, yerr=np.sqrt(data_hist), fmt="o", label="data", color="black")
    ax.plot(pos, VBF_hist, label="VBF", color="yellow", linestyle="--", linewidth=2)
    print("vbf purity", VBF_hist.sum()/sig_hist.sum())
    ax.plot(pos, sig_hist, label="signal", color="red", linewidth=3)
    
    ax.set_xlabel(r"$m_{ll\gamma}$ [GeV/c$^2$]")
    ax.set_ylabel("Events")
    ax.set_xlim(100, 180)
    ax.legend()
    
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/figs/hmass_dis_{i}.png")