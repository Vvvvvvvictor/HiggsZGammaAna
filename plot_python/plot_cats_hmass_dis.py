import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import os
plt.style.use(hep.style.CMS)

# 定义输出目录
output_dir = "/eos/home-j/jiehan/root/outputs_rui_run2/"

# 添加颜色字典
color_dict = {r"Z$+\gamma$": "#3f90da", "Z+Fake Photon": "#ffa90e", r"VBSZ+$\gamma$": "#92dadd", 
              r"t$\bar{t}$": "#e76300", r"t$\gamma$/t$\bar{t}\gamma$": "#bd1f01", 
              "multiboson": "#832db6", r"t$\bar{t}$+X": "#94a4a2"}

# 添加背景标签映射
background_labels = {
    "ZGToLLG": r"Z$+\gamma$", 
    "DYJetsToLL": "Z+Fake Photon", 
    "EWKZ2J": r"VBSZ+$\gamma$"
}

def read_config(file_name):
    with open(file_name, "r") as f:
        return list(map(float, f.readline().strip().split(" ")))
  
def get_hist(file_name, tree_name, var_name, bin_boundaries, ranges):
    data = uproot.open(file_name)[tree_name].arrays([var_name, "weight", "bdt_score_t", "n_jets"], library="pd").query(f"bdt_score_t>{bin_boundaries[0]} & bdt_score_t<{bin_boundaries[1]}")
    weight = data["weight"]
    data = data[var_name]
    return np.histogram(data, weights=weight, range=ranges[1:], bins=ranges[0])

def get_err_hist(file_name, tree_name, var_name, bin_boundaries, ranges):
    data = uproot.open(file_name)[tree_name].arrays([var_name, "weight", "bdt_score_t", "n_jets"], library="pd").query(f"bdt_score_t>{bin_boundaries[0]} & bdt_score_t<{bin_boundaries[1]}")
    weight = data["weight"]
    data = data[var_name]
    return np.histogram(data, weights=weight**2, range=ranges[1:], bins=ranges[0])

boundaries = read_config(f"{output_dir}test/significances/bin_boundaries_1D_two_jet.txt")
backgrounds = ["ZGToLLG", "DYJetsToLL", "EWKZ2J"]
signal = ["ggH_M125", "VBF_M125"] #, "ZH_M125", "WplusH_M125", "WminusH_M125", "ttH_M125"]
data = ["Data"]

if not os.path.exists("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/"):
    os.makedirs("/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/")

folder = "test"  # 可以改为 "test" 以切换文件夹

for i in range(len(boundaries)-1):
    bkg_hist = []
    for bkg in backgrounds:
        hist, bins = get_hist(f"{output_dir}{folder}/two_jet/{bkg}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
        bkg_hist.append(hist)
        
    sig_hist = np.zeros(len(bins)-1)
    VBF_hist, bins = get_hist(f"{output_dir}{folder}/two_jet/VBF_M125.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
    # sig_err_hist = np.zeros(len(bins)-1)
    for sig in signal:
        hist, bins = get_hist(f"{output_dir}{folder}/two_jet/{sig}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
        sig_hist += hist
        # sig_err_hist += get_err_hist(f"{output_dir}{folder}/two_jet/{sig}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])[0]
        
    data_hist = np.zeros(len(bins)-1)
    for dat in data:
        hist, bins = get_hist(f"{output_dir}{folder}/two_jet/{dat}.root", "two_jet", "H_mass", boundaries[i:], [80, 100, 180])
        data_hist += hist

    sf = np.sum(data_hist) / np.sum(bkg_hist)
    bkg_hist = [bkg * sf for bkg in bkg_hist]
    
    fig, ax = plt.subplots()
    # 使用颜色字典和标签映射
    labels = [background_labels[bkg] for bkg in backgrounds]
    colors = [color_dict[background_labels[bkg]] for bkg in backgrounds]
    hep.histplot(bkg_hist, bins, stack=True, label=labels, histtype="fill", color=colors, ax=ax)
    mask = (bins[:-1] + bins[1:]) / 2
    data_hist_masked = np.where((mask > 120) & (mask < 130), 0, data_hist)
    ax.errorbar(mask, data_hist_masked, yerr=np.sqrt(data_hist_masked), fmt="o", label="data", color="black")
    ax.plot((bins[:-1] + bins[1:]) / 2, VBF_hist, label="VBF", color="yellow", linestyle="--", linewidth=2)
    ax.plot((bins[:-1] + bins[1:]) / 2, sig_hist, label="signal", color="red", linewidth=3)
    ax.set_xlim(100, 180)
    ax.set_ylim(0, max(np.sum(bkg_hist, axis=0).max(), data_hist.max()) * 1.1)
    ax.set_xlabel("Higgs Mass [GeV]")
    ax.set_ylabel("Events")
    ax.title.set_text(f"Events in category {i} of VBF channel")
    # ax.title.set_text(f"All events in VBF channel")
    ax.legend()
    # plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/hmass_dis_all.png")
    plt.savefig(f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/two_jet/hmass_dis_{i}.png")
    plt.clf()