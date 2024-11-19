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

boundaries = read_config("/eos/home-j/jiehan/root/outputs/significances/bin_boundaries_1D_two_jet.txt")
backgrounds = ["ZGToLLG", "DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]

if not os.path.exists("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/"):
    os.makedirs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/")

for i in range(4):
    bkg_hist = []
    for bkg in backgrounds:
        hist, bins = get_hist(f"/eos/user/j/jiehan/root/outputs/two_jet/{bkg}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
        bkg_hist.append(hist)
    
    fig, ax = plt.subplots()
    hep.histplot(bkg_hist, bins, stack=True, label=backgrounds, histtype="fill", ax=ax)
    ax.set_xlim(100, 180)
    ax.set_xlabel("Higgs Mass [GeV]")
    ax.set_ylabel("Events")
    ax.legend()
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/hmass_dis_{i}.png")
    plt.clf()