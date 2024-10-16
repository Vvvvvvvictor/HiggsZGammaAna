import matplotlib.pyplot as plt
import uproot
import os
import numpy as np
import mplhep as hep
from pdb import set_trace

plt.style.use(hep.style.CMS)

# Settings
file_path = "/eos/home-j/jiehan/root"
type_dirs = ["", "_run2_gdr_4", "_run2_jet25"]
signal_samples = ["ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"]
sig_lables = ["ggH", "VBF", r"W$^{-}$H", r"W$^{+}$H", "ZH", "ttH"]
background_samples = ["ZGToLLG", "DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
bkg_lables = ["SM ZG", "DY", "EWK Z", "EWK ZG", "WG", "TT", "TTG", "TG", "WW", "WZ", "ZZ", "ttZ", "ttW"]
trees = ["zero_to_one_jet", "two_jet", "VH", "ZH", "ttH_lep", "ttH_had"]

# Load the root file
def load_yields(samples, file_path, type_dirs, trees, name):
    yields = []
    for type in type_dirs:
        for sample in samples:
            sample_yield = np.zeros(len(trees))
            if name == "signal":
                sample_path = f"{file_path}/skimmed_ntuples{type}/{sample}_M125"
            else:
                sample_path = f"{file_path}/skimmed_ntuples{type}/{sample}"
            for file in os.listdir(sample_path):
                if file.endswith(".root"):
                    temp_yields = np.array([])
                    print(f"Opening {sample_path}/{file}")
                    f = uproot.open(f"{sample_path}/{file}")
                    for chan in trees:
                        data = f[chan].arrays(["H_mass", "weight"], library="pd")
                        temp_yields = np.append(temp_yields, np.sum(data.query("H_mass > 100")["weight"]))
                    sample_yield += temp_yields
            yields.append(sample_yield)
    return np.array(yields).reshape(len(type_dirs), len(samples), len(trees))

# Plotting function
def plot_yields(yields, samples, types, file_name, name):
    channels = np.arange(len(trees)) * 2.5
    n_samples = len(samples)
    n_datasets = len(yields)
    width = 0.4
    gap_between_sets = 0.1
    
    # Normalize yields
    normalized_yields = yields / yields[0].sum(axis=0)
    print("normalized_yields", np.sum(normalized_yields, axis=1))

    # 定义20种颜色
    colors = ['blue', 'green', 'red', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'magenta', 'yellow', 'black', 'lime', 'teal', 'navy', 'maroon', 'gold', 'indigo', 'tan']
    
    fig, ax = plt.subplots()
    
    for i, dataset in enumerate(normalized_yields):
        x = channels + i * (width + gap_between_sets)
        bottom = np.zeros(len(channels))
        
        for j in range(n_samples):
            ax.bar(x, dataset[j], width, bottom=bottom, label=samples[j] if i == 0 else "", color=colors[j])
            bottom += dataset[j]

    for i in range(n_datasets):
        mid_x = channels[0] + i * (width + gap_between_sets)
        ypos = np.sum(normalized_yields[i], axis=0)[0]
        ax.text(mid_x, 1.05 * ypos, types[i], ha='center', va='bottom', fontsize=12)
        ax.annotate('', xy=(mid_x, 1.02 * ypos), xytext=(mid_x, 1.05 * ypos), arrowprops=dict(facecolor='black', shrink=0.05))

    ax.set_xticks(channels + (n_datasets - 1) * (width + gap_between_sets) / 2)
    ax.set_xticklabels(["01J", "VBF", "VH", "ZH", "ttH(lep)", "ttH(had)"])
    ax.set_ylim(0, np.max(np.sum(normalized_yields, axis=1)) * (1.2 if name == "signal" else 1.33))
    ax.legend(ncol=4, loc='upper right')
    ax.set_xlabel('Channels')
    ax.set_ylabel('Rel. Yields')
    ax.grid(axis='y')
    
    hep.cms.label("Preliminary", ax=ax)
    plt.tight_layout()
    plt.savefig(file_name)
    plt.clf()

# Define types for datasets
types = ["dR3", "dR4", "jet veto"]

# Load yields
sig_yields = load_yields(signal_samples, file_path, type_dirs, trees, "signal")
bkg_yields = load_yields(background_samples, file_path, type_dirs, trees, "background")

# Plot yields
plot_yields(sig_yields, sig_lables, types, "signal_yields_change_with_cond.png", "signal")
plot_yields(bkg_yields, bkg_lables, types, "background_yields_change_with_cond.png", "background")
