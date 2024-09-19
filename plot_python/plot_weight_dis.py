import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

def plot_weight(file_name):
    input_file = uproot.open(file_name)
    for tree in ["zero_jet", "one_jet", "two_jet", "ZH", "ttH_had", "ttH_lep"]:
        weight = input_file[tree].arrays(["weight"], library="pd")["weight"]
        weight_point = np.unique(weight)
        weight_dis = np.array([np.sum(weight[weight==point]) for point in weight_point])
        plt.scatter(weight_point, weight_dis, label=tree, s=3)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Weight")
        plt.ylabel("Number of events")
        plt.legend()
        plt.savefig("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/WminusH"+file_name.split("/")[-1].replace(".root", f"_{tree}.png"))
        plt.clf()

plot_weight("/eos/home-j/jiehan/root/skimmed_ntuples_run3/WminusH_M125/2022postEE.root")