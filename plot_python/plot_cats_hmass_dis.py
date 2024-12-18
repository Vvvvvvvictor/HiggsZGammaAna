import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import os
plt.style.use(hep.style.CMS)

def read_config(file_name):
    with open(file_name, "r") as f:
        return list(map(float, f.readline().strip().split(" ")))
  
def get_hist(file_name, tree_name, var_name, bin_boundaries, ranges):
    data = uproot.open(file_name)[tree_name].arrays([var_name, "weight", "bdt_score_t"], library="pd").query(f"bdt_score_t>{bin_boundaries[0]} & bdt_score_t<{bin_boundaries[1]}")
    weight = data["weight"]
    data = data[var_name]
    return np.histogram(data, weights=weight, range=ranges[1:], bins=ranges[0])

def get_err_hist(file_name, tree_name, var_name, bin_boundaries, ranges):
    data = uproot.open(file_name)[tree_name].arrays([var_name, "weight", "bdt_score_t"], library="pd").query(f"bdt_score_t>{bin_boundaries[0]} & bdt_score_t<{bin_boundaries[1]}")
    weight = data["weight"]
    data = data[var_name]
    return np.histogram(data, weights=weight**2, range=ranges[1:], bins=ranges[0])

boundaries = read_config("/eos/home-j/jiehan/root/outputs/significances/bin_boundaries_1D_two_jet.txt")
backgrounds = ["ZGToLLG", "DYJetsToLL", "ZG2JToG2L2J"]
signal = ["ggH_M125", "VBF_M125", "ZH_M125", "WplusH_M125", "WminusH_M125", "ttH_M125"]
data = ["Data"]

if not os.path.exists("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/"):
    os.makedirs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/")

for i in range(len(boundaries)-1):
    bkg_hist = []
    for bkg in backgrounds:
        hist, bins = get_hist(f"/eos/user/j/jiehan/root/outputs/two_jet/{bkg}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
        bkg_hist.append(hist)
        
    sig_hist = np.zeros(len(bins)-1)
    VBF_hist, bins = get_hist(f"/eos/user/j/jiehan/root/outputs/two_jet/VBF_M125.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
    # sig_err_hist = np.zeros(len(bins)-1)
    for sig in signal:
        hist, bins = get_hist(f"/eos/user/j/jiehan/root/outputs/two_jet/{sig}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
        sig_hist += hist
        # sig_err_hist += get_err_hist(f"/eos/user/j/jiehan/root/outputs/two_jet/{sig}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])[0]
        
    data_hist = np.zeros(len(bins)-1)
    for dat in data:
        hist, bins = get_hist(f"/eos/user/j/jiehan/root/outputs/two_jet/{dat}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
        data_hist += hist
    
    fig, ax = plt.subplots()
    hep.histplot(bkg_hist, bins, stack=True, label=backgrounds, histtype="fill", ax=ax)
    mask = (bins[:-1] + bins[1:]) / 2
    data_hist_masked = np.where((mask > 120) & (mask < 130), 0, data_hist)
    ax.errorbar(mask, data_hist_masked, yerr=np.sqrt(data_hist_masked), fmt="o", label="data", color="black")
    ax.plot((bins[:-1] + bins[1:]) / 2, VBF_hist, label="VBF", color="yellow", linestyle="--", linewidth=2)
    ax.plot((bins[:-1] + bins[1:]) / 2, sig_hist, label="signal", color="red", linewidth=3)
    ax.set_xlim(100, 180)
    ax.set_xlabel("Higgs Mass [GeV]")
    ax.set_ylabel("Events")
    ax.legend()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/hmass_dis_{i}.png")
    plt.clf()